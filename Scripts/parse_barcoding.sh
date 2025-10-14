#!/bin/bash

#SBATCH --partition=larsms
#SBATCH --mem=30G

R1=$1
R2=$2

#placeholder genome_dir bc i'm doing alignment on all files together later
split-pipe \
    --mode pre \
    --chemistry v1 \
    --genome_dir  ../../../scRNA-DNA/20231023_SplitTAP_pilot/homebrewParse/genomes/hg38 \
    --fq1 $R1 \
    --fq2 $R2 \
    --output_dir Data/barcoded_fastqs
