#!/bin/bash
#SBATCH --partition=larsms
#SBATCH --time=5:20:00
#SBATCH --cpus-per-task=20

ml load biology
ml load samtools/1.16.1
ml load star/2.7.10b

fastq=$1
#suffix=$(echo "$fastq" | grep -oP '(?<=_)\d+(?=\.fastq\.gz)')

#STAR --runMode genomeGenerate \
#    --genomeDir STAR/full_transcripts \
#     --genomeFastaFiles STAR/full_transcripts.fa \
#     --sjdbGTFfile STAR/full_transcripts.gtf \
#     --sjdbOverhang 149 \
#     --genomeSAindexNbases 5

# Align RNA-seq reads to the genome
STAR --runThreadN 20 \
     --genomeDir STAR/full_transcripts \
     --readFilesIn $fastq \
     --sjdbGTFfile STAR/full_transcripts.gtf \
     --outSAMmultNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix Data/BAM/ \
     --chimSegmentMin 20 \
     --outFilterMultimapNmax 30 \
     --chimOutType Junctions SeparateSAMold \
     --outFilterScoreMinOverLread .5 \
     --outFilterMatchNminOverLread .5 \
     --readFilesCommand zcat \
     --outReadsUnmapped Fastx

