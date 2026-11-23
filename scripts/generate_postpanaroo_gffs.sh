#!/bin/bash
#SBATCH --job-name=generate_postpanaroo_gffs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=panaroo_generate_gffs_%j.log
#SBATCH --error=panaroo_generate_gffs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=lab-colibri-hmem
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

set -euo pipefail
source "$(dirname "$0")/../config.sh"

INPUT_DIR="$POSTPANAROO_INPUT_DIR"
PANAROO_OUT="$PANAROO_RESULTS_DIR"
SCRIPT="$(dirname "$0")/generate_postpanaroo_gff_files.py"

echo "======================================="
echo "Starting custom Panaroo GFF generation"
echo "Date:   $(date)"
echo "Host:   $(hostname)"
echo "Input:  $INPUT_DIR"
echo "Output: $PANAROO_OUT"
echo "Script: $SCRIPT"
echo "======================================="

module load miniconda3

# Use the Panaroo environment
conda activate "$PANAROO_CONDA_ENV"

# Prevent contamination from ~/.local packages
export PYTHONNOUSERSITE=1

mkdir -p "$PANAROO_OUT"

echo "Python:"
which python

echo "Running GFF generation..."

python "$SCRIPT" \
    -o "$PANAROO_OUT" \
    -i "$INPUT_DIR"/*.gff \
    -f prokka \
    -t 1

echo "======================================="
echo "Finished: $(date)"
echo "======================================="
