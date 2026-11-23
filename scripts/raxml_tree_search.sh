#!/bin/bash
#SBATCH --job-name=raxml_search
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evho@ucsc.edu
#SBATCH --output=raxml_search_%j.log
#SBATCH --error=raxml_search_%j.err
#SBATCH --mem=200G
#SBATCH --cpus-per-task=56
#SBATCH --partition=lab-colibri-hmem
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module load raxml-ng/1.2.0
source "$(dirname "$0")/../config.sh"

OUT="$RAXML_RESULTS_DIR"

raxml-ng --search1 \
    --msa ${OUT}/core_genome_dedup_parse.raxml.rba \
    --model GTR+G \
    --prefix ${OUT}/core_genome_search \
    --threads 56 \
    --seed 12345
