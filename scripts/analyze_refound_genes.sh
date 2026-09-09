#!/bin/bash
source "$(dirname "$0")/../config.sh"

ALIGN_DIR="$PROJECT_ROOT/panaroo_results/ref_guided_DNA_alignments"
OUTPUT="$PROJECT_ROOT/panaroo_results/refound_length_report.txt"

echo -e "file\ttotal_seqs\tmodal_length\trefound_diff_length\trefound_diff_ids" > "$OUTPUT"

for f in "$ALIGN_DIR"/*.aln.fas; do
    basename=$(basename "$f")

    # Get total sequence count
    total=$(grep -c "^>" "$f")

    # Get modal (most common) length
    modal=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' "$f" | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')

    # Find refound sequences with different length
    diff_ids=$(awk -v modal="$modal" '/^>/{if(seq && /refound/ && length(seq) != modal) print header; header=$0; seq=""} !/^>/{seq=seq$0} END{if(seq && header ~ /refound/ && length(seq) != modal) print header}' "$f" | tr '\n' ',' | sed 's/,$//')

    # Count them
    if [ -z "$diff_ids" ]; then
        count=0
    else
        count=$(echo "$diff_ids" | tr ',' '\n' | wc -l)
    fi

    echo -e "$basename\t$total\t$modal\t$count\t$diff_ids" >> "$OUTPUT"
done

echo "Done! Report written to $OUTPUT"
