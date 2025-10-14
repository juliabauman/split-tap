#!/bin/bash

# Usage check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.fastq.gz output.fastq.gz"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
LINKER="LINKER"

zcat "$INPUT" | awk -v linker="$LINKER" '
NR%4 == 1 {header=$0; next}
NR%4 == 2 {
    linker_idx = index($0, linker)
    if (linker_idx > 0) {
        spacer = substr($0, linker_idx + length(linker))
        seq = substr($0, 1, linker_idx - 1)
        header = header ",SPACER1=" spacer
        spacer_len = length(spacer) + length(linker)
    } else {
        header = header ",SPACER1="
        seq = $0
        spacer_len = 0
    }
    print header
    print seq
    next
}
NR%4 == 3 {print; next}
NR%4 == 0 {
    if (spacer_len > 0) {
        print substr($0, 1, length($0) - spacer_len)
    } else {
        print $0
    }
}
' | gzip > "$OUTPUT"

echo "✅ Processed: $INPUT → $OUTPUT"

