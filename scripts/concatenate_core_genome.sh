#!/bin/bash
#SBATCH --job-name=core_genome_concat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=concatenate_core_genomes_%j.log
#SBATCH --error=concatenate_core_genomes_%j.err
#SBATCH --mem=16G

# Activate environment
module load miniconda3
source "$(dirname "$0")/../config.sh"
conda activate "$ALIGNMENT_ENV"

PANAROO_DIR="$PANAROO_RESULTS_DIR"

python "$(dirname "$0")/concatenate_core_genome.py" \
    -o ${PANAROO_DIR} \
    --filtered-core ${PANAROO_DIR}/ref_guided_DNA_alignments \
    --core_threshold 0.95 \
    --prefix ref_guided_core_gene
