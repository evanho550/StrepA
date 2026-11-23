import os
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq

from panaroo.isvalid import is_valid_folder


def get_options():
    import argparse

    description = 'Generate GFF annotation files for each isolate from the pan-genome graph'
    parser = argparse.ArgumentParser(description=description,
                                     prog='generate_gff_files')

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory",
                         type=lambda x: is_valid_folder(parser, x))

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files used to run panaroo",
        type=str,
        nargs='+')

    core = parser.add_argument_group('GFFs')
    core.add_argument(
        "-f",
        "--f",
        dest="format",
        choices={'prokka', 'gff3'},
        default='prokka')

    parser.add_argument("-t", "--threads", dest="n_cpu", type=int, default=1)
    parser.add_argument("--verbose", dest="verbose", action='store_true')

    return parser.parse_args()


def parse_all_gffs(list_of_isolate_names, input_list, verbose):
    ordered_parsed_gffs = []

    if len(input_list) == 1:
        with open(input_list[0], 'r') as inhandle:
            input_list = inhandle.read().splitlines()

    if verbose:
        list_of_isolate_names = tqdm(list_of_isolate_names)

    for isolate_id in list_of_isolate_names:
        parsed_gff = {}

        base_input = [os.path.basename(x).replace(".gff", "") for x in input_list]
        input_file = next((s[1] for s in zip(base_input, input_list) if isolate_id == s[0]))

        with open(input_file, 'r') as raw_file:
            lines = raw_file.read()   # ❗ FIX: removed .replace(',', '')

        split = lines.split('##FASTA')

        parsed_gff["fasta"] = split[1] if len(split) > 1 else None

        header, body = [], []
        for line in split[0].splitlines():
            if "##" in line:
                header.append(line)
            else:
                body.append(line)

        parsed_gff["header"] = "\n".join(header)
        parsed_gff["body"] = body

        ordered_parsed_gffs.append(parsed_gff)

    return ordered_parsed_gffs


def parse_gff_body(gff_body_list):
    parsed_gff = []

    for gff_line in gff_body_list:
        cols = gff_line.split("\t")
        if len(cols) < 9:
            continue

        if cols[2] != "CDS":
            continue

        parsed = {
            "seqid": cols[0],
            "type": cols[2],
            "start": cols[3],
            "end": cols[4],
            "score": cols[5],
            "strand": cols[6],
            "phase": cols[7],
        }

        attributes = ";".join(cols[8:]).split(";")
        for a in attributes:
            if a.startswith("ID="):
                parsed["ID"] = a.split("=")[1]
            elif "eC_number=" in a:
                parsed["eC_number"] = a.split("=")[1]
            elif "inference=" in a:
                parsed["inference"] = a.split("=")[1]
            elif "locus_tag=" in a:
                parsed["Locus_tag"] = a.split("=")[1]

        parsed_gff.append(parsed)

    return parsed_gff


def process_refound_gene(refound_id, pangenome_id, parsed_gff, refound_seqs,
                         output_dir, G):

    refound_seq = refound_seqs[refound_id]

    scaffold_id = refound_seq[0]
    start = refound_seq[1]
    stop = refound_seq[2]
    strand = refound_seq[3]

    gene_name = G.nodes[pangenome_id]["annotation"]
    gene_description = G.nodes[pangenome_id]["description"]
    has_paralog = "True" if G.nodes[pangenome_id]["paralog"] == 1 else "False"

    gff_attributes = ";".join([
        "ID=" + refound_id,
        "locus_tag=" + refound_id,
        "name=" + gene_name,
        "description=" + gene_description,
        "inference=panaroo refound gene",
        "has_pangenome_paralog=" + has_paralog
    ])

    return "\t".join([
        scaffold_id,
        "Panaroo_refound",
        "candidate_gene",
        str(start + 1),
        str(stop),
        ".",
        strand,
        ".",
        gff_attributes
    ])


def create_new_gffs(isolate_index, parsed_gffs, pp_isolate_genes,
                    gene_name_dic, refound_seqs, outdir, gff_format, G):

    new_gff_body_lines = []
    new_gff_body_lines.append(parsed_gffs[isolate_index]["header"])

    parsed_original_gffbody = parse_gff_body(parsed_gffs[isolate_index]["body"])

    if str(isolate_index) not in pp_isolate_genes:
        return []

    for pangenome_gene in pp_isolate_genes[str(isolate_index)]:

        for gene in pp_isolate_genes[str(isolate_index)][pangenome_gene]:

            if "refound" in gene:
                new_gff_body_lines.append(
                    process_refound_gene(
                        gene, pangenome_gene,
                        parsed_gffs[isolate_index],
                        refound_seqs, outdir, G
                    )
                )
                continue

            expected_id = gene_name_dic[gene]

            # --------------------------------------------------
            # First try exact ID matching (original Panaroo logic)
            # --------------------------------------------------
            original_gene_data = [
                x for x in parsed_original_gffbody
                if x.get("ID") == expected_id
            ]

            # --------------------------------------------------
            # Special-case fallback for programmed frameshift
            # genes stored as:
            #
            #     cds-XXXX
            #     cds-XXXX_1
            #
            # in gene_data.csv but only:
            #
            #     cds-XXXX
            #
            # in the original GFF.
            # --------------------------------------------------
            if len(original_gene_data) == 0 and expected_id.endswith("_1"):

                fallback_id = expected_id.rsplit("_", 1)[0]

                original_gene_data = [
                    x for x in parsed_original_gffbody
                    if x.get("ID") == fallback_id
                ]

                if len(original_gene_data) > 0:
                    print(
                        f"Frameshift fallback used: "
                        f"{expected_id} -> {fallback_id}"
                    )

            # --------------------------------------------------
            # Still missing after fallback
            # --------------------------------------------------
            if len(original_gene_data) == 0:

                print("\n========== MISSING GENE ==========")
                print("isolate_index:", isolate_index)
                print("panaroo gene:", gene)
                print("expected ID :", expected_id)

                if expected_id.endswith("_1"):
                    print("fallback ID :", expected_id.rsplit("_", 1)[0])

                print("==================================\n")

                raise RuntimeError(
                    f"Could not find {expected_id} in isolate {isolate_index}"
                )

            # --------------------------------------------------
            # Multiple matches after fallback.
            # This appears to occur for programmed
            # frameshift genes represented by multiple
            # CDS segments in the original GFF.
            # --------------------------------------------------
            if len(original_gene_data) > 1:

                print(
                    f"Frameshift CDS detected: "
                    f"{expected_id} ({len(original_gene_data)} segments)"
                )

                for hit in original_gene_data:
                    print(
                        "  ",
                        hit.get("ID"),
                        hit.get("seqid"),
                        hit.get("start"),
                        hit.get("end")
                    )

                # Use the leftmost CDS segment
                original_gene_data = sorted(
                    original_gene_data,
                    key=lambda x: int(x["start"])
                )

            # Convert list of matches into a single CDS record
            original_gene_data = original_gene_data[0]

            gene_name = G.nodes[pangenome_gene]["annotation"] or "No_name"
            gene_description = G.nodes[pangenome_gene]["description"]
            has_paralog = "True" if G.nodes[pangenome_gene]["paralog"] == 1 else "False"

            new_attributes = ";".join([
                "ID=" + expected_id,
                "name=" + gene_name,
                "locus_tag=" + original_gene_data.get("Locus_tag", "Panaroo_refound"),
                "description=" + gene_description,
                "pangenome_id=" + str(pangenome_gene),
                "panaroo_ID=" + gene,
                "eC_number=" + original_gene_data.get("eC_number", "None"),
                "prepanaroo_inference=" + original_gene_data.get("inference", "NA"),
                "has_pangenome_paralog=" + has_paralog
            ])

            new_gff_body_lines.append("\t".join([
                original_gene_data["seqid"],
                "Panaroo",
                original_gene_data["type"],
                original_gene_data["start"],
                original_gene_data["end"],
                original_gene_data["score"],
                original_gene_data["strand"],
                original_gene_data["phase"],
                new_attributes
            ]))

    new_gff_body_lines[1:] = sorted(
        new_gff_body_lines[1:],
        key=lambda x: (x.split('\t')[0], int(x.split('\t')[3]))
    )

    if gff_format == "prokka":
        new_gff_body_lines.append("##FASTA")
        new_gff_body_lines.append(parsed_gffs[isolate_index]["fasta"].strip())

    return new_gff_body_lines


def output_gff(isolate_name, gff_lines, outdir):
    with open(outdir + "postpanaroo_gffs/" + isolate_name + "_panaroo.gff", "w") as f:
        f.write("\n".join(gff_lines))


def main():
    args = get_options()
    args.output_dir = os.path.join(args.output_dir, "")

    if args.verbose:
        print("Loading panaroo output and input gff files...")

    seen = set()
    isolate_names = []
    gene_names = {}
    refound_seqs = {}

    with open(args.output_dir + "gene_data.csv") as infile:
        next(infile)
        for line in infile:
            s = line.strip().split(",")
            iso = s[0]
            gene_names[s[2]] = s[3]

            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)

            if "refound" in s[2]:
                loc, strand = s[-1].split(';')
                loc = loc.split(':')[1].split('-')
                strand = strand.split(':')[1]
                refound_seqs[s[2]] = (s[1], int(loc[0]), int(loc[1]), strand)

    G = nx.read_gml(args.output_dir + "final_graph.gml")

    parsed_gffs = parse_all_gffs(isolate_names, args.input_files, args.verbose)

    isolate_genes = {}

    for node in tqdm(G.nodes):
        for iso in G.nodes[node]["genomeIDs"].split(";"):
            isolate_genes.setdefault(iso, {})
            seqIDs = G.nodes[node]["seqIDs"]

            if isinstance(seqIDs, list):
                for sid in seqIDs:
                    if sid.split("_")[0] == iso:
                        isolate_genes[iso].setdefault(node, []).append(sid)
            else:
                if seqIDs.split("_")[0] == iso:
                    isolate_genes[iso].setdefault(node, []).append(seqIDs)

    if args.verbose:
        print("Creating new gff files...")

    new_gffs = Parallel(n_jobs=args.n_cpu, prefer="threads")(
        delayed(create_new_gffs)(
            i, parsed_gffs, isolate_genes,
            gene_names, refound_seqs,
            args.output_dir, args.format, G
        )
        for i in tqdm(range(len(isolate_names)))
    )

    os.makedirs(args.output_dir + "postpanaroo_gffs", exist_ok=True)

    for i, gff in enumerate(new_gffs):
        output_gff(isolate_names[i], gff, args.output_dir)


if __name__ == "__main__":
    main()
