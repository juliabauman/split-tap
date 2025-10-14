#!/bin/bash
#SBATCH --job-name=ftp_wget
#SBATCH --output=ftp_wget_%j.log
#SBATCH --time=6:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

wget --continue -r -nH --cut-dirs=1 \
  --user=Julia.Bauman \
  --password=nQ9WKz \
  ftp://38.122.175.98:2223/25095-01-03252025_191850

