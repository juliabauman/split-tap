#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --partition=larsms

ml load biology
ml load py-pandas

#for sample in $(ls *_L007_R1_001.fastq.gz | sed 's/_L007_R1_001.fastq.gz//'); do
#  echo "Merging R1 for $sample"
#  cat ${sample}_L007_R1_001.fastq.gz ${sample}_L008_R1_001.fastq.gz > ${sample}_L009_R1_001.fastq.gz
#
#  echo "Merging R2 for $sample"
#  cat ${sample}_L007_R2_001.fastq.gz ${sample}_L008_R2_001.fastq.gz > ${sample}_L009_R2_001.fastq.gz
#done

#rm *_L007_*.fastq.gz *_L008_*.fastq.gz

python3 ../Scripts/rename_fastqs.py

#make a backup tar file of all fastqs
tar -cvf /oak/stanford/groups/larsms/Users/jrbauman/backup_data/JB_ST_002/screen_all_samples.tar.gz -C . $(find . -maxdepth 1 -type d -regex './[A-C][0-9]+_[12]')
