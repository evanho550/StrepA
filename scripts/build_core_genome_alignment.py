#!/usr/bin/env python3
"""
build_core_genome_alignment.py

Builds two synteny-ordered concatenated core genome alignments for Gubbins:
  1. ALL core genes:       Corekaburra segments + GML-placed singleton genes
  2. COREKABURRA ONLY:     Genes confirmed in syntenic blocks only

Gene ordering logic:
  - Within segments:   consensus synteny order from Corekaburra core_segments.tsv
  - Between segments:  Panaroo final_graph.gml edge weights (genome support count)
  - Singleton genes:   placed by GML adjacency to nearest segment endpoint
  - Circular genome linearised by breaking at the minimum-weight inter-segment edge

Outputs (all written to --out_dir):
  {prefix}_all_alignment.fasta           Gubbins input — all core genes
  {prefix}_all_alignment.embl            Gubbins annotation header
  {prefix}_all_gene_order.tsv            Annotated gene order list
  {prefix}_corekaburra_alignment.fasta   Gubbins input — synteny-confirmed genes only
  {prefix}_corekaburra_alignment.embl
  {prefix}_corekaburra_gene_order.tsv
  {prefix}_diagnostic_report.txt        Full diagnostic report

Usage:
    python3 build_core_genome_alignment.py \
        --panaroo_dir     /path/to/panaroo_results \
        --corekaburra_dir /path/to/corekaburra_results_2 \
        --aln_dir         /path/to/ref_guided_DNA_alignments \
        --out_dir         /path/to/output \
        [--core_threshold 0.95] \
        [--aln_ext        .aln.fas] \
        [--prefix         strep_core]

─────────────────────────────────────────────────────────────────────────────
PSEUDOCODE (for notes)
─────────────────────────────────────────────────────────────────────────────

LOAD isolates from gene_data.csv  (N ≈ 1892)
LOAD Panaroo graph G from final_graph.gml

IDENTIFY core genes
  → nodes where node.size / N >= threshold (default 0.95)
  → note genes with hasEnd=1 (at contig edges in ≥1 genome; useful for diagnostics)

LOAD Corekaburra segments from core_segments.tsv
  → segment_to_genes: {segment → [gene1, gene2, ..., geneK]}  (Corekaburra position order)
  → gene_to_segment:  {gene → segment}

LOAD per-gene alignments from ref_guided_DNA_alignments/
  → for each core gene, open {gene}.aln.fas
  → warn and skip genes with no file
  → if duplicate isolate IDs, keep the sequence with fewer gaps

BUILD SUPER-NODES  (atomic units for graph traversal)
  For each Corekaburra segment:
    → one super-node; genes locked in Corekaburra order
    → endpoints = (first_gene, last_gene)  ← the two genes that face outward
  For each core gene NOT in any segment  (singleton):
    → one super-node containing just that one gene

BUILD SEGMENT-LEVEL GRAPH
  For each GML edge where BOTH endpoints are core genes in DIFFERENT super-nodes:
    → add edge between those super-nodes
    → weight = GML edge 'size'  (number of genomes supporting that adjacency)
    → if multiple GML edges connect the same pair, keep the highest-weight one
    → record which specific gene pair forms the connection  (needed for orientation)
  Internal segment edges (both genes in the same super-node) → skip; already ordered

FIND TRAVERSAL ORDER  (linearise the circular chromosome)
  Find minimum-weight edge in segment graph → break the circle here
  Start greedy traversal from one side of the break:
    → at each node, move to the highest-weight unvisited neighbor
    → determine each segment's orientation from the connecting gene:
         connecting gene == segment.first_gene  →  traverse forward
         connecting gene == segment.last_gene   →  traverse reverse
    → log ambiguous points (multiple unvisited neighbors with similar weights)
  Any super-nodes not reached  →  "unplaced"  (appended at end with warning)

EXPAND TRAVERSAL TO GENE LIST
  For each (super-node, orientation) in traversal:
    → emit genes in forward or reversed Corekaburra order
  Append unplaced super-nodes at end (labelled)

  Produce two lists:
    gene_order_all          all core genes in chromosome order
    gene_order_corekaburra  filter gene_order_all to corekaburra_segment genes only




FOR EACH GENE ORDER  (all / corekaburra-only):
  Build per-isolate concatenated alignment:
    → for each gene in order: append alignment sequence
    → if isolate missing a gene: append gap string of correct length
  Write:
    → {prefix}_{label}_alignment.fasta   (Gubbins FASTA input)
    → {prefix}_{label}_alignment.embl    (Gubbins annotation header)
    → {prefix}_{label}_gene_order.tsv    (position, gene, source, supernode, length)

WRITE DIAGNOSTIC REPORT
  → run parameters
  → summary counts (isolates, core genes, segments, singletons, missing, unplaced)
  → cycle break edge and its genome support weight
  → missing alignment files (with hasEnd flag)
  → Corekaburra segments table  (name, size, gene list)
  → singleton placement details  (segment-graph neighbours, max/min edge weights, hasEnd)
  → unplaced super-nodes
  → ambiguous traversal points
  → full ordered gene lists with source annotation for both outputs
─────────────────────────────────────────────────────────────────────────────
"""

import os
import sys
import csv
import argparse
import networkx as nx
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_isolates(gene_data_csv):
    """Return ordered list of unique isolate IDs from Panaroo gene_data.csv."""
    seen, isolates = set(), []
    with open(gene_data_csv) as fh:
        next(fh)  # skip header
        for line in fh:
            iso = line.split(",")[0].strip()
            if iso and iso not in seen:
                isolates.append(iso)
                seen.add(iso)
    return isolates


def get_core_genes(G, threshold, num_isolates):
    """
    Identify core gene nodes from the Panaroo GML graph.

    Returns:
        node_to_name:  {node_id: gene_name}  for nodes at or above threshold
        has_end_genes: set of gene names where hasEnd=1 (at contig edge in >=1 genome)
    """
    node_to_name = {}
    has_end_genes = set()
    for node in G.nodes():
        size = float(G.nodes[node].get("size", 0))
        if size / num_isolates >= threshold:
            name = G.nodes[node].get("name", str(node))
            node_to_name[node] = name
            if int(G.nodes[node].get("hasEnd", 0)) == 1:
                has_end_genes.add(name)
    return node_to_name, has_end_genes


def load_corekaburra_segments(core_segments_tsv):
    """
    Parse Corekaburra core_segments.tsv.

    Returns:
        segment_to_genes: {segment_name: [gene1, gene2, ...]}  in position order
        gene_to_segment:  {gene_name: segment_name}
    """
    raw = defaultdict(dict)
    with open(core_segments_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            seg  = row['Segment_name'].strip()
            pos  = int(row['Segment_position'])
            gene = row['Core_gene'].strip()
            raw[seg][pos] = gene

    segment_to_genes, gene_to_segment = {}, {}
    for seg, pos_dict in raw.items():
        ordered = [pos_dict[p] for p in sorted(pos_dict)]
        segment_to_genes[seg] = ordered
        for gene in ordered:
            gene_to_segment[gene] = seg
    return segment_to_genes, gene_to_segment


def load_alignments(aln_dir, core_gene_names, aln_ext):
    """
    Load per-gene FASTA alignments for all core genes.

    Returns:
        alignments:   {gene_name: {isolate_id: sequence_str}}
        gene_lengths: {gene_name: alignment_length_bp}
        missing_aln:  set of gene names with no alignment file found
    """
    alignments, gene_lengths, missing_aln = {}, {}, set()
    for gene in core_gene_names:
        path = os.path.join(aln_dir, gene + aln_ext)
        if not os.path.exists(path):
            missing_aln.add(gene)
            continue
        records = list(SeqIO.parse(path, 'fasta'))
        if not records:
            missing_aln.add(gene)
            continue
        gene_dict = {}
        gene_len  = len(records[0].seq)
        for rec in records:
            iso = rec.id.split(";")[0].strip()
            seq = str(rec.seq)
            # Keep sequence with fewer gaps if same isolate appears twice
            if iso in gene_dict:
                if seq.count("-") < gene_dict[iso].count("-"):
                    gene_dict[iso] = seq
            else:
                gene_dict[iso] = seq
        alignments[gene]   = gene_dict
        gene_lengths[gene] = gene_len
    return alignments, gene_lengths, missing_aln


# ─────────────────────────────────────────────────────────────────────────────
# 2. SUPER-NODE CONSTRUCTION
# ─────────────────────────────────────────────────────────────────────────────

def build_supernodes(segment_to_genes, core_gene_names, gene_to_segment):
    """
    Build super-nodes — the atomic units the segment-graph traversal works with.

    Each Corekaburra segment becomes one super-node (genes locked in Corekaburra order).
    Each core gene not in any segment becomes a singleton super-node.

    Returns:
        supernodes:        {sn_id: {'type', 'genes', 'endpoints', 'name'}}
        gene_to_supernode: {gene_name: sn_id}
    """
    supernodes, gene_to_supernode = {}, {}

    # ── Corekaburra segment super-nodes ──────────────────────────────────────
    for seg_name, genes in segment_to_genes.items():
        core_g = [g for g in genes if g in core_gene_names]
        if not core_g:
            continue
        sn_id = f"SEG:{seg_name}"
        supernodes[sn_id] = {
            'type':      'segment',
            'genes':     core_g,
            'endpoints': (core_g[0], core_g[-1]),  # genes that face outward
            'name':      seg_name,
        }
        for gene in core_g:
            gene_to_supernode[gene] = sn_id

    # ── Singleton super-nodes ─────────────────────────────────────────────────
    for gene in core_gene_names:
        if gene not in gene_to_segment:
            sn_id = f"SING:{gene}"
            supernodes[sn_id] = {
                'type':      'singleton',
                'genes':     [gene],
                'endpoints': (gene, gene),
                'name':      gene,
            }
            gene_to_supernode[gene] = sn_id

    return supernodes, gene_to_supernode


# ─────────────────────────────────────────────────────────────────────────────
# 3. SEGMENT-LEVEL GRAPH
# ─────────────────────────────────────────────────────────────────────────────

def build_segment_graph(G, node_to_name, supernodes, gene_to_supernode):
    """
    Build a weighted graph between super-nodes using Panaroo GML adjacency edges.

    For each GML edge where both endpoints are core genes in DIFFERENT super-nodes:
      - Add/update an edge between those super-nodes
      - Weight = GML edge 'size' (number of genomes supporting that adjacency)
      - Keep only the highest-weight edge between any pair of super-nodes
      - Record the connecting gene pair (used later for orientation)

    Internal segment edges (both genes in the same super-node) are skipped;
    their order is already determined by Corekaburra.

    Returns:
        SG:              networkx.Graph with super-node edges
        edge_gene_pairs: {(sn_a, sn_b) sorted tuple: (gene_a, gene_b, weight)}
    """
    SG = nx.Graph()
    SG.add_nodes_from(supernodes.keys())
    edge_gene_pairs = {}

    for u, v, data in G.edges(data=True):
        # Skip if either endpoint is not a core gene
        if u not in node_to_name or v not in node_to_name:
            continue
        gene_u = node_to_name[u]
        gene_v = node_to_name[v]

        # Skip if gene not mapped to a super-node (shouldn't happen, but guard)
        if gene_u not in gene_to_supernode or gene_v not in gene_to_supernode:
            continue

        sn_u = gene_to_supernode[gene_u]
        sn_v = gene_to_supernode[gene_v]

        # Skip internal edges — order already fixed by Corekaburra
        if sn_u == sn_v:
            continue

        weight   = int(data.get('size', 1))
        edge_key = tuple(sorted([sn_u, sn_v]))

        if SG.has_edge(sn_u, sn_v):
            if weight > SG[sn_u][sn_v]['weight']:
                SG[sn_u][sn_v]['weight'] = weight
                edge_gene_pairs[edge_key] = (gene_u, gene_v, weight)
        else:
            SG.add_edge(sn_u, sn_v, weight=weight)
            edge_gene_pairs[edge_key] = (gene_u, gene_v, weight)

    return SG, edge_gene_pairs


# ─────────────────────────────────────────────────────────────────────────────
# 4. TRAVERSAL
# ─────────────────────────────────────────────────────────────────────────────

def _get_orientation(prev_sn, curr_sn, supernodes, edge_gene_pairs):
    """
    Determine whether curr_sn should be traversed forward or reverse,
    based on which endpoint gene connects to prev_sn.

      connecting gene == segment.endpoints[0]  →  forward  (enter from front)
      connecting gene == segment.endpoints[1]  →  reverse  (enter from back)
      connecting gene is internal              →  forward  (fallback, with warning)
    """
    if prev_sn is None:
        return 'forward'
    edge_key = tuple(sorted([prev_sn, curr_sn]))
    if edge_key not in edge_gene_pairs:
        return 'forward'

    gene_a, gene_b, _ = edge_gene_pairs[edge_key]
    sn_data = supernodes[curr_sn]

    # Find which of the two genes belongs to curr_sn
    if gene_a in sn_data['genes']:
        connecting = gene_a
    elif gene_b in sn_data['genes']:
        connecting = gene_b
    else:
        return 'forward'  # shouldn't happen

    ep_first, ep_last = sn_data['endpoints']
    if connecting == ep_first:
        return 'forward'
    elif connecting == ep_last:
        return 'reverse'
    else:
        # Connecting gene is internal to segment — rearrangement artifact
        # Default to forward; flagged in diagnostic report
        print(f"    [WARNING] Internal gene '{connecting}' connects "
              f"'{prev_sn}' to '{curr_sn}' — defaulting to forward orientation")
        return 'forward'


def find_traversal_order(SG, supernodes, edge_gene_pairs):
    if SG.number_of_nodes() == 0:
        return [], [], [], None

    # Global minimum-weight edge (for diagnostic report only)
    break_edge = None
    min_w = float('inf')
    for u, v, d in SG.edges(data=True):
        w = d.get('weight', 0)
        if w < min_w:
            min_w = w
            break_edge = (u, v)

    visited, traversal, ambiguous_log = set(), [], []
    components = sorted(nx.connected_components(SG), key=len, reverse=True)

    for component in components:
        # Minimum-weight edge within this component → cycle-break start point
        comp_break = None
        comp_min_w = float('inf')
        for u, v, d in SG.edges(data=True):
            if u in component and v in component:
                w = d.get('weight', 0)
                if w < comp_min_w:
                    comp_min_w = w
                    comp_break = (u, v)

        current = comp_break[0] if comp_break else next(iter(component))
        prev    = None

        # Keep restarting within this component until every node is visited
        while True:
            while current is not None:
                visited.add(current)
                orientation = _get_orientation(prev, current, supernodes, edge_gene_pairs)
                traversal.append((current, orientation))

                candidates = [
                    (nb, SG[current][nb].get('weight', 0))
                    for nb in SG.neighbors(current)
                    if nb not in visited
                ]
                candidates.sort(key=lambda x: -x[1])

                if not candidates:
                    current = None
                else:
                    if len(candidates) > 1:
                        ambiguous_log.append({
                            'supernode':     current,
                            'chosen':        candidates[0][0],
                            'chosen_weight': candidates[0][1],
                            'alternatives':  candidates[1:],
                        })
                    prev, current = current, candidates[0][0]

            # Dead end — check if any nodes in this component are still unvisited
            remaining = [n for n in component if n not in visited]
            if not remaining:
                break          # whole component done
            prev    = None     # fresh start: no orientation inheritance
            current = remaining[0]

    unplaced = [sn for sn in SG.nodes() if sn not in visited]
    return traversal, ambiguous_log, unplaced, break_edge


def expand_traversal(traversal, supernodes, unplaced):
    gene_order = []
    seen_genes = set()

    for i, (sn_id, orientation) in enumerate(traversal):
        sn    = supernodes[sn_id]
        genes = sn['genes'] if orientation == 'forward' else list(reversed(sn['genes']))
        src   = 'corekaburra_segment' if sn['type'] == 'segment' else 'singleton_gml'

        start = 0
        if i > 0 and genes:
            prev_sn, prev_orient = traversal[i - 1]
            prev_genes = (
                supernodes[prev_sn]['genes']
                if prev_orient == 'forward'
                else list(reversed(supernodes[prev_sn]['genes']))
            )
            if prev_genes and genes[0] == prev_genes[-1]:
                start = 1

        for gene in genes[start:]:
            if gene not in seen_genes:
                gene_order.append((gene, src, sn['name']))
                seen_genes.add(gene)

    for sn_id in unplaced:
        sn  = supernodes[sn_id]
        src = 'unplaced_segment' if sn['type'] == 'segment' else 'unplaced_singleton'
        for gene in sn['genes']:
            if gene not in seen_genes:
                gene_order.append((gene, src, sn['name']))
                seen_genes.add(gene)

    return gene_order


# ─────────────────────────────────────────────────────────────────────────────
# 5. ALIGNMENT CONCATENATION
# ─────────────────────────────────────────────────────────────────────────────

def concatenate_alignments(gene_order, alignments, gene_lengths, isolates):
    """
    Concatenate per-gene alignments in gene_order for all isolates.
    Genes with no alignment file are silently skipped (already warned at load time).
    Isolates missing a gene are gap-padded to the correct alignment length.

    Returns:
        concat_seqs:    {isolate: concatenated_sequence_str}
        used_gene_order: [(gene, source, sn_name)] genes actually included
    """
    used = [(g, s, n) for g, s, n in gene_order if g in alignments]
    concat_seqs = {iso: "" for iso in isolates}

    for gene, _, _ in used:
        gene_len = gene_lengths[gene]
        aln      = alignments[gene]
        for iso in isolates:
            concat_seqs[iso] += aln.get(iso, "-" * gene_len)

    return concat_seqs, used


# ─────────────────────────────────────────────────────────────────────────────
# 6. OUTPUT WRITING
# ─────────────────────────────────────────────────────────────────────────────

def write_fasta(concat_seqs, isolates, output_path):
    """Write concatenated alignment as FASTA."""
    records = [
        SeqRecord(Seq(concat_seqs[iso]), id=iso, description="")
        for iso in isolates
    ]
    SeqIO.write(records, output_path, "fasta")
    print(f"      Written: {output_path}")


def write_embl_header(used_gene_order, gene_lengths, output_path):
    """Write EMBL-format annotation header for Gubbins."""
    total_len = sum(gene_lengths[g] for g, _, _ in used_gene_order if g in gene_lengths)
    pos = 1
    with open(output_path, 'w') as fh:
        fh.write(f"ID   Genome standard; DNA; PRO; {total_len} BP.\nXX\n")
        fh.write("FH   Key             Location/Qualifiers\nFH\n")
        for gene, _, _ in used_gene_order:
            if gene not in gene_lengths:
                continue
            end = pos + gene_lengths[gene] - 1
            fh.write(
                f"FT   feature         {pos}..{end}\n"
                f"FT                   /label={gene}\n"
                f"FT                   /locus_tag={gene}\n"
            )
            pos += gene_lengths[gene]
        fh.write("XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T; 0 other;\n//\n")
    print(f"      Written: {output_path}")


def write_gene_order_tsv(used_gene_order, gene_lengths, output_path):
    """Write annotated gene order as TSV for reference and manual inspection."""
    with open(output_path, 'w') as fh:
        fh.write("position\tgene\tsource\tsupernode\tlength_bp\n")
        for i, (gene, source, sn_name) in enumerate(used_gene_order, 1):
            gene_len = gene_lengths.get(gene, 'NA')
            fh.write(f"{i}\t{gene}\t{source}\t{sn_name}\t{gene_len}\n")
    print(f"      Written: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. DIAGNOSTIC REPORT
# ─────────────────────────────────────────────────────────────────────────────

def write_diagnostic_report(output_path, args,
                             num_isolates, num_core_genes,
                             segment_to_genes, supernodes,
                             gene_order_all, gene_order_ck,
                             alignments, gene_lengths, missing_aln,
                             has_end_genes, ambiguous_log, unplaced,
                             SG, edge_gene_pairs, break_edge):

    n_seg_genes  = sum(len(v) for v in segment_to_genes.values())
    n_singletons = sum(1 for sn in supernodes.values() if sn['type'] == 'singleton')
    n_all        = sum(1 for g, _, _ in gene_order_all if g in alignments)
    n_ck         = sum(1 for g, _, _ in gene_order_ck  if g in alignments)

    with open(output_path, 'w') as fh:

        def w(line=""):
            fh.write(line + "\n")

        w("=" * 72)
        w("CORE GENOME ALIGNMENT — DIAGNOSTIC REPORT")
        w("=" * 72)

        # ── Run parameters ─────────────────────────────────────────────────
        w("\n── RUN PARAMETERS ───────────────────────────────────────────────────")
        w(f"  Panaroo directory:           {args.panaroo_dir}")
        w(f"  Corekaburra directory:       {args.corekaburra_dir}")
        w(f"  Alignment directory:         {args.aln_dir}")
        w(f"  Core threshold:              {args.core_threshold:.0%}")
        w(f"  Alignment extension:         {args.aln_ext}")
        w(f"  Output prefix:               {args.prefix}")

        # ── Summary ────────────────────────────────────────────────────────
        w("\n── SUMMARY ──────────────────────────────────────────────────────────")
        w(f"  Isolates:                              {num_isolates}")
        w(f"  Core genes (>={args.core_threshold:.0%}):                {num_core_genes}")
        w(f"  Corekaburra segments:                  {len(segment_to_genes)}")
        w(f"  Genes in Corekaburra segments:         {n_seg_genes}")
        w(f"  Singleton genes (GML-placed):          {n_singletons}")
        w(f"  Genes with hasEnd=1 (contig-edge):     {len(has_end_genes)}")
        w(f"  Missing alignment files:               {len(missing_aln)}")
        w(f"  Unplaced super-nodes:                  {len(unplaced)}")
        w(f"  Ambiguous traversal points:            {len(ambiguous_log)}")
        w()
        w(f"  ALL alignment  — genes included:       {n_all}")
        w(f"  COREKABURRA alignment — genes:         {n_ck}")

        if break_edge:
            bw = SG[break_edge[0]][break_edge[1]].get('weight', '?')
            be_names = (
                supernodes[break_edge[0]]['name'],
                supernodes[break_edge[1]]['name'],
            )
            w()
            w(f"  Cycle break between:  {be_names[0]}  <-->  {be_names[1]}")
            w(f"  Break edge weight:    {bw} / {num_isolates} genomes")

        # ── Missing alignments ─────────────────────────────────────────────
        if missing_aln:
            w("\n── MISSING ALIGNMENT FILES ──────────────────────────────────────────")
            w("  Core genes excluded — no .aln.fas file found:")
            for g in sorted(missing_aln):
                tag = "  [hasEnd=1]" if g in has_end_genes else ""
                w(f"    {g}{tag}")

        # ── Corekaburra segments ───────────────────────────────────────────
        w("\n── COREKABURRA SEGMENTS (largest first) ─────────────────────────────")
        w(f"  {'Segment':<35}  {'Size':>5}  Genes")
        w("  " + "-" * 70)
        for seg, genes in sorted(segment_to_genes.items(), key=lambda x: -len(x[1])):
            w(f"  {seg:<35}  {len(genes):>5}  {', '.join(genes)}")

        # ── Singleton placement details ─────────────────────────────────────
        w("\n── SINGLETON GENES — GML PLACEMENT DETAILS ──────────────────────────")
        w(f"  {'Gene':<25}  {'Neighbours':>10}  {'Max wt':>8}  {'Min wt':>8}  hasEnd")
        w("  " + "-" * 65)
        for sn_id, sn in sorted(supernodes.items()):
            if sn['type'] != 'singleton':
                continue
            gene = sn['genes'][0]
            nbs  = list(SG.neighbors(sn_id))
            he   = "yes" if gene in has_end_genes else "no"
            if nbs:
                wts = [SG[sn_id][nb].get('weight', 0) for nb in nbs]
                w(f"  {gene:<25}  {len(nbs):>10}  {max(wts):>8}  {min(wts):>8}  {he}")
            else:
                w(f"  {gene:<25}  {'0':>10}  {'N/A':>8}  {'N/A':>8}  {he}")

        # ── Unplaced super-nodes ───────────────────────────────────────────
        if unplaced:
            w("\n── UNPLACED SUPER-NODES (disconnected from main chromosome) ─────────")
            w("  Appended at end of ALL alignment; excluded from COREKABURRA alignment.")
            for sn_id in unplaced:
                sn = supernodes[sn_id]
                w(f"  [{sn['type']}]  {sn['name']}:  {', '.join(sn['genes'])}")

        # ── Ambiguous traversal points ─────────────────────────────────────
        if ambiguous_log:
            w("\n── AMBIGUOUS TRAVERSAL POINTS ───────────────────────────────────────")
            w("  Points where multiple routes existed (highest weight was chosen):")
            for entry in ambiguous_log:
                w(f"  At:     {entry['supernode']}")
                w(f"    Chose:  {entry['chosen']}  (weight {entry['chosen_weight']})")
                for alt, wt in entry['alternatives']:
                    w(f"    Alt:    {alt}  (weight {wt})")

        # ── Full ordered gene list — ALL ───────────────────────────────────
        w("\n── GENE ORDER: ALL CORE GENES ───────────────────────────────────────")
        w(f"  {'Pos':>5}  {'Gene':<25}  {'Source':<22}  Supernode")
        w("  " + "-" * 80)
        pos = 1
        for gene, source, sn_name in gene_order_all:
            if gene in alignments:
                he = " *" if gene in has_end_genes else ""
                w(f"  {pos:>5}  {gene:<25}  {source:<22}  {sn_name}{he}")
                pos += 1
        w("  (* = gene flagged hasEnd=1 in Panaroo graph)")

        # ── Full ordered gene list — COREKABURRA ONLY ─────────────────────
        w("\n── GENE ORDER: COREKABURRA ONLY ─────────────────────────────────────")
        w(f"  {'Pos':>5}  {'Gene':<25}  Segment")
        w("  " + "-" * 60)
        pos = 1
        for gene, source, sn_name in gene_order_ck:
            if gene in alignments:
                w(f"  {pos:>5}  {gene:<25}  {sn_name}")
                pos += 1

        w("\n" + "=" * 72)
        w("END OF REPORT")

    print(f"      Written: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 8. MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build synteny-ordered core genome alignments for Gubbins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--panaroo_dir",      required=True,
                        help="Panaroo output dir (contains final_graph.gml and gene_data.csv)")
    parser.add_argument("--corekaburra_dir",  required=True,
                        help="Corekaburra output dir (contains core_segments.tsv)")
    parser.add_argument("--aln_dir",          required=True,
                        help="Directory of per-gene FASTA alignments")
    parser.add_argument("--out_dir",          required=True,
                        help="Output directory (created if it does not exist)")
    parser.add_argument("--core_threshold",   type=float, default=0.95,
                        help="Core gene presence threshold")
    parser.add_argument("--aln_ext",          default=".aln.fas",
                        help="Alignment file extension")
    parser.add_argument("--prefix",           default="strep_core",
                        help="Output file prefix")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Step 1: Load isolates ───────────────────────────────────────────────
    print("[1/8] Loading isolates from gene_data.csv...")
    isolates = load_isolates(os.path.join(args.panaroo_dir, "gene_data.csv"))
    print(f"      {len(isolates)} isolates found")

    # ── Step 2: Load Panaroo graph ──────────────────────────────────────────
    print("[2/8] Loading Panaroo graph (this may take a moment)...")
    G = nx.read_gml(os.path.join(args.panaroo_dir, "final_graph.gml"))
    print(f"      {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ── Step 3: Identify core genes ─────────────────────────────────────────
    print("[3/8] Identifying core genes...")
    node_to_name, has_end_genes = get_core_genes(G, args.core_threshold, len(isolates))
    core_gene_names = set(node_to_name.values())
    print(f"      {len(core_gene_names)} core genes at {args.core_threshold:.0%} threshold")
    print(f"      {len(has_end_genes)} genes with hasEnd=1 (at contig edges in >=1 genome)")

    # ── Step 4: Load Corekaburra segments ──────────────────────────────────
    print("[4/8] Loading Corekaburra segments...")
    segment_to_genes, gene_to_segment = load_corekaburra_segments(
        os.path.join(args.corekaburra_dir, "core_segments.tsv")
    )
    n_seg_genes = sum(len(v) for v in segment_to_genes.values())
    print(f"      {len(segment_to_genes)} segments covering {n_seg_genes} genes")

    # ── Step 5: Load gene alignments ───────────────────────────────────────
    print("[5/8] Loading gene alignments...")
    alignments, gene_lengths, missing_aln = load_alignments(
        args.aln_dir, core_gene_names, args.aln_ext
    )
    print(f"      {len(alignments)} alignments loaded")
    if missing_aln:
        print(f"      [WARNING] {len(missing_aln)} core genes have no alignment file — will be excluded")

    # ── Step 6: Build super-nodes and segment graph ─────────────────────────
    print("[6/8] Building super-nodes and segment-level graph...")
    supernodes, gene_to_supernode = build_supernodes(
        segment_to_genes, core_gene_names, gene_to_segment
    )
    n_seg_sn  = sum(1 for sn in supernodes.values() if sn['type'] == 'segment')
    n_sing_sn = sum(1 for sn in supernodes.values() if sn['type'] == 'singleton')
    print(f"      {n_seg_sn} segment super-nodes, {n_sing_sn} singleton super-nodes")

    SG, edge_gene_pairs = build_segment_graph(G, node_to_name, supernodes, gene_to_supernode)
    print(f"      Segment graph: {SG.number_of_nodes()} nodes, {SG.number_of_edges()} edges")

    # ── Step 7: Traverse and expand ─────────────────────────────────────────
    print("[7/8] Traversing segment graph to determine gene order...")
    traversal, ambiguous_log, unplaced, break_edge = find_traversal_order(
        SG, supernodes, edge_gene_pairs
    )
    print(f"      {len(traversal)} super-nodes placed, {len(unplaced)} unplaced")
    if ambiguous_log:
        print(f"      {len(ambiguous_log)} ambiguous traversal points (see diagnostic report)")
    if unplaced:
        print(f"      [WARNING] {len(unplaced)} super-node(s) disconnected — appended at end of ALL alignment")

    gene_order_all = expand_traversal(traversal, supernodes, unplaced)

    # Corekaburra-only: keep segment genes only
    gene_order_ck = [
        (g, s, n) for g, s, n in gene_order_all
        if s == 'corekaburra_segment'
    ]

    # ── Step 8: Concatenate and write ───────────────────────────────────────
    print("[8/8] Concatenating alignments and writing outputs...")

    for label, gene_order, suffix in [
        ("ALL",         gene_order_all, "all"),
        ("COREKABURRA", gene_order_ck,  "corekaburra"),
    ]:
        print(f"\n  [{label}]")
        concat_seqs, used = concatenate_alignments(
            gene_order, alignments, gene_lengths, isolates
        )
        total_bp = sum(gene_lengths[g] for g, _, _ in used)
        print(f"    {len(used)} genes  |  {total_bp:,} bp total")

        # Per-isolate gap-padding warnings
        n_warned = 0
        for iso in isolates:
            n_miss = sum(1 for g, _, _ in used if iso not in alignments.get(g, {}))
            if n_miss and n_warned < 10:
                pct = 100 * n_miss / len(used)
                print(f"    [WARNING] {iso}: missing {n_miss}/{len(used)} ({pct:.1f}%) genes — gap-padded")
                n_warned += 1
        if n_warned == 10:
            print(f"    [WARNING] ... (further gap-padding warnings suppressed; see diagnostic report)")

        base = os.path.join(args.out_dir, f"{args.prefix}_{suffix}")
        write_fasta(concat_seqs, isolates, base + "_alignment.fasta")
        write_embl_header(used, gene_lengths, base + "_alignment.embl")
        write_gene_order_tsv(used, gene_lengths, base + "_gene_order.tsv")

    # ── Diagnostic report ────────────────────────────────────────────────────
    print("\n  [REPORT]")
    report_path = os.path.join(args.out_dir, f"{args.prefix}_diagnostic_report.txt")
    write_diagnostic_report(
        report_path, args,
        len(isolates), len(core_gene_names),
        segment_to_genes, supernodes,
        gene_order_all, gene_order_ck,
        alignments, gene_lengths, missing_aln,
        has_end_genes, ambiguous_log, unplaced,
        SG, edge_gene_pairs, break_edge,
    )

    print(f"\n[DONE] All outputs written to: {args.out_dir}")


if __name__ == "__main__":
    main()

