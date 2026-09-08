#!/bin/bash
#SBATCH --job-name=ref_align
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=panaroo_results/logs/aln_%A_%a.out
#SBATCH --error=panaroo_results/logs/aln_%A_%a.err
#SBATCH --array=1-63          # Fill in after generating file_list.txt
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source "$(dirname "$0")/../config.sh"
cd "$PROJECT_ROOT/panaroo_results"

# Safety check: Make sure task ID is within file list bounds
num_files=$(wc -l < failed_genes_file_list.txt)
if [ "$SLURM_ARRAY_TASK_ID" -gt "$num_files" ]; then
  echo "Task ID $SLURM_ARRAY_TASK_ID exceeds number of files ($num_files). Exiting."
  exit 1
fi

# Read the file corresponding to the current array task ID (1-based)
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" failed_genes_file_list.txt)

# Extract the base name without directory and extension
basename=$(basename "$input_file" .fasta)

# Look up protein alignment — check both singleton (.aa.fasta) and mafft (.mafft.aa.fasta)
COMBINED_DIR="$PROJECT_ROOT/panaroo_results/combined_aligned_protein_sequences"
if [ -f "${COMBINED_DIR}/${basename}.mafft.aa.fasta" ]; then
    aa_aligned="${COMBINED_DIR}/${basename}.mafft.aa.fasta"
elif [ -f "${COMBINED_DIR}/${basename}.aa.fasta" ]; then
    aa_aligned="${COMBINED_DIR}/${basename}.aa.fasta"
else
    echo "ERROR: No protein alignment found for ${basename}. Exiting."
    exit 1
fi

# Define output path
dna_output="$PROJECT_ROOT/panaroo_results/ref_guided_DNA_alignments/${basename}.aln.fas"

module load miniconda3
conda activate "$ALIGNMENT_ENV"

# Run the alignment script
python "$(dirname "$0")/ref_guided_alignment.py" \    
	"$input_file" "$aa_aligned" "$dna_output"
