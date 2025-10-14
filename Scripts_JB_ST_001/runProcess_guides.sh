#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=8 
#SBATCH --partition=larsms

ml load biology
ml load python/3.9
ml load py-numpy/1.24.2_py39
ml load py-pandas/2.0.1_py39
ml load py-biopython/1.79_py39
ml load py-regex/20247.24_py39
ml load py-schwimmbad/0.3.2_py39

#python3 Scripts/process_guides_new.py Data/barcoded_fastqs/all_CROP_barcoded.fastq.gz Data/guides.txt

python3 Scripts/theseus_impl.py Data/guides.txt Data/guides_filtered.txt 0.2

python3 Scripts/count_guides.py Data/guides_filtered.txt Data/guide_counts.txt
