#!/bin/bash
module load miniconda3
source "$(dirname "$0")/../config.sh"
conda activate "$ALIGNMENT_ENV"

ALIGN_DIR="$PROJECT_ROOT/panaroo_results/ref_guided_DNA_alignments"
REPORT="$PROJECT_ROOT/panaroo_results/refound_length_report.txt"

while IFS=$'\t' read -r filename total modal count ids; do
    [ "$filename" = "file" ] && continue
    [ "$count" -eq 0 ] && continue

    filepath="$ALIGN_DIR/$filename"

    python3 - "$filepath" "$ids" <<'EOF'
import sys
from Bio import SeqIO

input_file = sys.argv[1]
bad_ids = set(sys.argv[2].split(','))
# Strip leading > from IDs if present
bad_ids = {i.lstrip('>').split()[0] for i in bad_ids}

records = [r for r in SeqIO.parse(input_file, "fasta") if r.id not in bad_ids]

with open(input_file, 'w') as out:
    SeqIO.write(records, out, "fasta")

print(f"Modified: {input_file} (removed {len(bad_ids)} refound sequences)")
EOF

done < "$REPORT"

echo "Done!"
