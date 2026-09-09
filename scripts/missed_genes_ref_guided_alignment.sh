#!/bin/bash
module load miniconda3
source "$(dirname "$0")/../config.sh"
conda activate "$ALIGNMENT_ENV"

COMBINED_DIR="$PROJECT_ROOT/panaroo_results/combined_aligned_protein_sequences"
DNA_ALIGN_DIR="$PROJECT_ROOT/panaroo_results/ref_guided_DNA_alignments"
SCRIPT="$(dirname "$0")/ref_guided_alignment.py"
FILE_LIST="$PROJECT_ROOT/panaroo_results/all_gene_sequences_file_list.txt"

for gene in gloA group_212 group_343 group_510 group_75; do
    dna=$(grep "/${gene}\.fasta" "$FILE_LIST")
    protein=$(ls "${COMBINED_DIR}/${gene}.aa.fasta" \
                 "${COMBINED_DIR}/${gene}.mafft.aa.fasta" 2>/dev/null | head -1)
    output="${DNA_ALIGN_DIR}/${gene}.aln.fas"

    if [ -z "$dna" ]; then
        echo "ERROR: DNA file not found for $gene"
        continue
    fi
    if [ -z "$protein" ]; then
        echo "ERROR: Protein file not found for $gene"
        continue
    fi

    echo "Running: $gene"
    echo "  DNA:     $dna"
    echo "  Protein: $protein"
    echo "  Output:  $output"
    python "$SCRIPT" "$dna" "$protein" "$output"
done

echo "Done!"
