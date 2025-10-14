#!/bin/bash
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8

ml load system
ml load pigz/2.4

INPUT=$1
OUTPUT=$2
INDEX=$3

zcat "$INPUT" | awk -v i7="$INDEX" 'NR % 4 == 1 {
    split($0, parts, "__")
    if (length(parts) >= 5) {
        parts[4] = parts[4] "_" i7
        $0 = parts[1]
        for (i = 2; i <= 5; i++) $0 = $0 "__" parts[i]
    }
}
{ print }' | pigz -p 8 > "$OUTPUT"

