#!/bin/bash

#SBATCH --partition=larsms

STAR --runMode genomeGenerate \
    --genomeDir STAR/full_transcripts \
     --genomeFastaFiles STAR/full_transcripts.fa \
     --sjdbGTFfile STAR/full_transcripts.gtf \
     --sjdbOverhang 149 \
     --genomeSAindexNbases 5
