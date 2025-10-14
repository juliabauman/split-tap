#!/bin/bash

#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8
#SBATCH --mem=28G
#SBATCH --time=3:00:00

ml load biology
ml load samtools

#samtools cat -o Data/featureCounts/all_sample/mid_merge.bam \
#  Data/featureCounts/A10_body_chunk_aa.bam.featureCounts.bam \
#  Data/featureCounts/A10_body_chunk_ab.bam.featureCounts.bam \
#  Data/featureCounts/A10_body_chunk_ac.bam.featureCounts.bam \
#  Data/featureCounts/A10_body_chunk_ad.bam.featureCounts.bam \
#  Data/featureCounts/A10_body_chunk_ae.bam.featureCounts.bam \
#  Data/featureCounts/A11_body_chunk_aa.bam.featureCounts.bam \
#  Data/featureCounts/A11_body_chunk_ab.bam.featureCounts.bam \
#  Data/featureCounts/A11_body_chunk_ac.bam.featureCounts.bam \
#  Data/featureCounts/A11_body_chunk_ad.bam.featureCounts.bam \
#  Data/featureCounts/A11_body_chunk_ae.bam.featureCounts.bam \
#  Data/featureCounts/A12_body_chunk_aa.bam.featureCounts.bam \
#  Data/featureCounts/A12_body_chunk_ab.bam.featureCounts.bam \
#  Data/featureCounts/A12_body_chunk_ac.bam.featureCounts.bam \
#  Data/featureCounts/A12_body_chunk_ad.bam.featureCounts.bam \
#  Data/featureCounts/A12_body_chunk_ae.bam.featureCounts.bam \
#  Data/featureCounts/B1_body_chunk_aa.bam.featureCounts.bam \
#  Data/featureCounts/B1_body_chunk_ab.bam.featureCounts.bam \
#  Data/featureCounts/B1_body_chunk_ac.bam.featureCounts.bam \
#  Data/featureCounts/B1_body_chunk_ad.bam.featureCounts.bam \
#  Data/featureCounts/B1_body_chunk_ae.bam.featureCounts.bam


#Try different sort parameters
samtools sort -n -@ 8 -o Data/featureCounts/all_sample/mid_merge_sorted.bam Data/featureCounts/all_sample/mid_merge.bam
samtools index Data/featureCounts/all_sample/mid_merge_sorted.bam

# Before sort
samtools view Data/featureCounts/all_sample/mid_merge.bam | \
  awk '{for(i=12;i<=NF;i++) if($i ~ /^XT:Z:/){sub(/^XT:Z:/,"",$i); print $i; break}}' | \
  sort | uniq -c | sort -nr > xt_tags_before.txt

# After sort
samtools view Data/featureCounts/all_sample/mid_merge_sorted.bam
  awk '{for(i=12;i<=NF;i++) if($i ~ /^XT:Z:/){sub(/^XT:Z:/,"",$i); print $i; break}}' | \
  sort | uniq -c | sort -nr > xt_tags_after.txt



