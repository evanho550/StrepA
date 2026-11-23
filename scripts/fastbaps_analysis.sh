#!/bin/bash
#SBATCH --job-name=fastbaps
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=fastbaps_%j.log
#SBATCH --error=fastbaps_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --partition=lab-colibri-hmem
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module load miniconda3
source "$(dirname "$0")/../config.sh"
conda activate "$ALIGNMENT_ENV"

Rscript "$(dirname "$0")/fastbaps_analysis.R"
