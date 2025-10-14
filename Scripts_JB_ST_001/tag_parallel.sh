#!/bin/bash

#SBATCH --array=1-20
#SBATCH --partition=larsms

ml load biology
ml load system
ml load ncurses/6.0
ml load py-biopython/1.79_py39
ml load py-pysam/0.18.0_py39
ml load samtools/1.16.1

input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cdna_input_bams.txt)
output=${input/.out.bam/_tagged.bam}

samtools index "$input"
python3 Scripts/add_tags.py "$input" "$output"

