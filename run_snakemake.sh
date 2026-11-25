#!/bin/bash

#SBATCH --partition=larsms
#SBATCH --time=48:00:00
#SBATCH --mem=100G

ml R

source /home/users/jrbauman/miniconda3/etc/profile.d/conda.sh
conda activate split-tap

snakemake --cores 8 --rerun-incomplete
