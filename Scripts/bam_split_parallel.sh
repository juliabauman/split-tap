#!/bin/bash
#SBATCH --partition=larsms
#SBATCH --time=6:00:00

# Directory containing the BAM files
INPUT_DIR="Data/star_outputs"
OUTPUT_DIR="Data/featureCounts"
LOG_DIR="logs"

# Ensure the output and log directories exist
#mkdir -p "$OUTPUT_DIR"
#mkdir -p "$LOG_DIR"

# Loop over each BAM file in the input directory
for bam in "$INPUT_DIR"/*/Aligned.sortedByCoord.out_tagged.bam; do
  base=$(basename "$bam" .bam)
  sample=$(basename $(dirname "$bam"))
  mkdir -p "${INPUT_DIR}/${sample}/chunked"
  sbatch --export=ALL,bam="$bam",base="$base",OUTPUT_DIR="$OUTPUT_DIR",LOG_DIR="$LOG_DIR" Scripts/ft_cts_parallel.sh
done

