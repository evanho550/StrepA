#!/bin/bash

# Load shared paths (edit config.sh with your real paths — see config.sh.example)
source "$(dirname "$0")/../config.sh"

mkdir -p "$REMOVED_GFF_DIR"

python3 - "$METADATA_FILE" <<'EOF' | tr -d '\r' > /tmp/checkm_lookup.tsv
import csv
import sys

metadata = sys.argv[1]

with open(metadata, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        accession = row["Assembly Accession"].strip()
        paired = row["Assembly Paired Assembly Accession"].strip() or "NA"
        checkm = row["CheckM completeness"].strip()
        print(f"{accession}\t{paired}\t{checkm}")
EOF

declare -A checkm_score
declare -A paired_acc

while IFS=$'\t' read -r accession paired checkm; do
    checkm_score["$accession"]="$checkm"
    paired_acc["$accession"]="$paired"
done < /tmp/checkm_lookup.tsv

moved=0
kept=0

for gff in "$GFF_DIR"/*.gff; do
    [[ -f "$gff" ]] || continue
    filename=$(basename "$gff")
    accession="${filename%.gff}"
    reason=""

    if [[ "$accession" == GCA_* ]]; then
        paired="${paired_acc[$accession]}"
        if [[ -n "$paired" && "$paired" != "NA" && -f "$GFF_DIR/${paired}.gff" ]]; then
            reason="GCA with existing paired GCF present ($paired)"
        fi
    fi

    if [[ -z "$reason" ]]; then
        score="${checkm_score[$accession]}"
        if [[ -z "$score" ]]; then
            reason="Missing or no CheckM completeness data"
        else
            passes=$(awk -v s="$score" 'BEGIN { print (s+0 >= 95) ? "yes" : "no" }')
            if [[ "$passes" == "no" ]]; then
                reason="CheckM completeness $score < 95"
            fi
        fi
    fi

    if [[ -n "$reason" ]]; then
        echo "MOVING:  $filename  ($reason)"
        mv "$gff" "$REMOVED_GFF_DIR/"
        ((moved++))
    else
        echo "KEEPING: $filename  (CheckM: ${checkm_score[$accession]})"
        ((kept++))
    fi
done

echo ""
echo "============================="
echo "Done.  Moved: $moved | Kept: $kept"
echo "============================="
