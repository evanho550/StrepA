#!/bin/bash
#SBATCH --job-name=panaroo_05_18_2026    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=evho@ucsc.edu   # Where to send mail
#SBATCH --output=panaroo_results/panaroo.out
#SBATCH --error=panaroo_results/panaroo.err
#SBATCH --mem=40G                    # Allocated 20 G memory
#SBATCH --nodes=1
module load panaroo
source "$(dirname "$0")/../config.sh"
cd "$PROJECT_ROOT/converted_gffs"
mash_database="$MASH_DB"
panaroo -i *.gff -o "$PANAROO_RESULTS_DIR" --clean-mode strict -f 0.9 -a pan --aligner none -t 1 --remove-invalid-genes
