import os
import argparse
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_core_gene_nodes(G, threshold, num_isolates):
    core_nodes = []
    for node in G.nodes():
        if float(G.nodes[node]["size"]) / float(num_isolates) >= threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(alignments_dir, core_names, output_dir, output_prefix):
    all_files = os.listdir(alignments_dir)

    core_filenames = sorted([
        x for x in all_files
        if x.endswith(".aln.fas") and x.replace(".aln.fas", "") in core_names
    ])

    aligned_names = set(x.replace(".aln.fas", "") for x in core_filenames)
    missing = set(core_names) - aligned_names
    if missing:
        print(f"[WARNING] {len(missing)} core genes have no alignment file and will be skipped:")
        for m in sorted(missing):
            print(f"  {m}")

    gene_alignments = []
    isolates = set()

    for filename in core_filenames:
        gene_name = filename.replace(".aln.fas", "")
        file_path = os.path.join(alignments_dir, filename)
        records = list(SeqIO.parse(file_path, 'fasta'))

        if not records:
            print(f"[WARNING] Empty file: {filename} — skipping.")
            continue

        gene_length = len(records[0].seq)

        gene_dict = {}
        for record in records:
            genome_id = record.id.split(";")[0]
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(gene_dict[genome_id]).count("-"):
                    gene_dict[genome_id] = record.seq
            else:
                gene_dict[genome_id] = record.seq
            isolates.add(genome_id)

        gene_alignments.append((gene_name, gene_dict, gene_length))

    if not gene_alignments:
        raise RuntimeError("No valid alignments found. Check your input directory and file names.")

    total_len = sum(g[2] for g in gene_alignments)
    print(f"\n[INFO] {len(gene_alignments)} core genes | {len(isolates)} isolates | {total_len} bp total")

    for iso in sorted(isolates):
        n_missing = sum(1 for g in gene_alignments if iso not in g[1])
        if n_missing:
            pct = 100 * n_missing / len(gene_alignments)
            print(f"[WARNING] {iso} missing {n_missing}/{len(gene_alignments)} ({pct:.1f}%) genes — gap-padded.")

    isolate_aln = []
    for iso in sorted(isolates):
        concat_seq = ""
        for gene_name, gene_dict, gene_length in gene_alignments:
            concat_seq += str(gene_dict[iso]) if iso in gene_dict else "-" * gene_length
        isolate_aln.append(SeqRecord(Seq(concat_seq), id=iso, description=""))

    aln_out = os.path.join(output_dir, output_prefix + "_alignment.aln")
    SeqIO.write(isolate_aln, aln_out, "fasta")
    print(f"[INFO] Alignment written to: {aln_out}")

    write_alignment_header(gene_alignments, output_dir, output_prefix)


def write_alignment_header(gene_alignments, output_dir, prefix):
    entries = []
    pos = 1
    for gene_name, gene_dict, gene_length in gene_alignments:
        end = pos + gene_length - 1
        entries.append(
            f"FT   feature         {pos}..{end}\n"
            f"FT                   /label={gene_name}\n"
            f"FT                   /locus_tag={gene_name}\n"
        )
        pos += gene_length

    total_len = pos - 1
    embl_out = os.path.join(output_dir, prefix + "_alignment_header.embl")
    with open(embl_out, "w") as fh:
        fh.write(f"ID   Genome standard; DNA; PRO; {total_len} BP.\nXX\n")
        fh.write("FH   Key             Location/Qualifiers\nFH\n")
        for entry in entries:
            fh.write(entry)
        fh.write("XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T; 0 other;\n//\n")
    print(f"[INFO] EMBL header written to: {embl_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate codon-aligned core gene FASTAs into a core genome MSA"
    )
    parser.add_argument("-o", "--out_dir", dest="output_dir", required=True,
                        help="Panaroo output dir (contains final_graph.gml and gene_data.csv)")
    parser.add_argument("--core_threshold", dest="core", type=float, default=0.95,
                        help="Core genome presence threshold (default: 0.95)")
    parser.add_argument("--filtered-core", dest="prefiltered", default=None,
                        help="Directory of codon-aligned .aln.fas files")
    parser.add_argument("--prefix", dest="prefix", default="core_genome",
                        help="Output file prefix (default: core_genome)")
    args = parser.parse_args()

    output_dir = args.output_dir.rstrip("/") + "/"

    # Always load isolate count from gene_data.csv
    seen = set()
    isolate_names = []
    with open(os.path.join(output_dir, "gene_data.csv")) as fh:
        next(fh)
        for line in fh:
            iso = line.split(",")[0]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)
    print(f"[INFO] {len(isolate_names)} isolates found in gene_data.csv")

    # Always load graph — needed for core filtering regardless of mode
    G = nx.read_gml(os.path.join(output_dir, "final_graph.gml"))
    core_nodes = get_core_gene_nodes(G, args.core, len(isolate_names))
    core_names = set(G.nodes[x]["name"] for x in core_nodes)
    print(f"[INFO] {len(core_names)} core genes at threshold {args.core}")

    if args.prefiltered:
        available = set(
            x.replace(".aln.fas", "")
            for x in os.listdir(args.prefiltered)
            if x.endswith(".aln.fas")
        )
        to_use = core_names & available
        not_aligned = core_names - available
        not_core = available - core_names

        print(f"[INFO] {len(to_use)} core genes have alignments and will be concatenated")
        if not_aligned:
            print(f"[WARNING] {len(not_aligned)} core genes have no alignment file and will be skipped:")
            for g in sorted(not_aligned):
                print(f"  {g}")
        if not_core:
            print(f"[INFO] {len(not_core)} aligned genes excluded (below core threshold)")

        concatenate_core_genome_alignments(args.prefiltered, to_use, output_dir, args.prefix)

    else:
        concatenate_core_genome_alignments(
            os.path.join(output_dir, "aligned_gene_sequences/"),
            core_names, output_dir, args.prefix
        )
