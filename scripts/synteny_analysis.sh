#!/bin/bash
#SBATCH --job-name=corekaburra
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=corekaburra_%j.log
#SBATCH --error=corekaburra_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=lab-colibri-hmem
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module load miniconda3
source "$(dirname "$0")/../config.sh"
conda activate "$COREKABURRA_ENV"
export PYTHONNOUSERSITE=1

Corekaburra \
    -ig "$POSTPANAROO_GFFS_DIR"/*.gff \
    -ip "$PANAROO_RESULTS_DIR"/ \
    -cc 0.95 \
    -o "$COREKABURRA_RESULTS_DIR" \
    -c 4
