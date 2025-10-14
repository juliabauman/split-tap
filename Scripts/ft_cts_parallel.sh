#!/bin/bash
#SBATCH --array=1-26
#SBATCH --partition=larsms

ml load biology
ml load samtools/1.16.1

input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" aligned_tagged_bams.txt)

bash Scripts/feat_cts.sh "$input"

