#!/bin/bash
#SBATCH --job-name=umi_spacer_collapse
#SBATCH --array=0-57
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --partition=larsms
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

ml load python/3.9
ml load py-pandas

mapfile -t CHUNK_FILES < <(ls Data/CROP_reads_chunk/*)
INPUT_FILE=${CHUNK_FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$INPUT_FILE" .tsv)
OUTPUT_FILE="Data/umi_collapse_temp/output_${BASENAME}.tsv"
UMI_LOG="Data/umi_collapse_temp/umi_corrections_${BASENAME}.tsv"
SPACER_LOG="Data/umi_collapse_temp/spacer_corrections_${BASENAME}.tsv"

python3 Scripts/collapse_umis.py $INPUT_FILE $OUTPUT_FILE $UMI_LOG $SPACER_LOG

