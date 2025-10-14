#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --partition=larsms
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

ml load biology
ml load samtools
ml load py-pandas



## === NOTE: All files' combined BAM is too big for featureCounts to work properly. Can't run the below. === ##


#INPUT=$1
#OUTPUT_DIR=$2

#featureCounts -a STAR/full_transcripts.gtf -o Data/amplicon_assigned -g gene_name -R BAM -M --primary -T 16 $INPUT

#mv "Data/Aligned_sorted_final.out.bam.featureCounts.bam" $OUTPUT_DIR

#samtools sort -@ 8 -m 6G "${OUTPUT_DIR}/Aligned_sorted_final.out.bam.featureCounts.bam" -o "${OUTPUT_DIR}/Aligned_featureCounts_sorted.bam"

#samtools index "${OUTPUT_DIR}/Aligned_featureCounts_sorted.bam"





### ==== trying featureCounts on the individual STAR aligned BAMs, split into even smaller chunks  ==== ###


bam=$1

base=$(basename "$(dirname "$bam")")

mkdir -p Data/star_outputs/${base}/chunked


echo "🔪 Splitting $bam into chunks..."
samtools view "$bam" | split -l 50000000 - "Data/star_outputs/${base}/chunked/${base}_body_chunk_"


for file in Data/star_outputs/${base}/chunked/${base}_body_chunk_*; do
  cat <(samtools view -H "$bam") "$file" | \
    samtools view -bS -o "${file}.bam" -
  rm "$file"
  echo "✅ Created chunk: ${file}.bam"
done


## Do featureCounts for each chunk
for chunk in Data/star_outputs/${base}/chunked/*.bam; do
  region=$(basename "$chunk" .bam)
  tag_base="${base}_${region}"
  bam_base=$(basename "$chunk")

  featureCounts -a STAR/full_transcripts.gtf -o Data/featureCounts/amplicon_assigned_${tag_base}.txt -g gene_name -R BAM -M --primary -T 8 $chunk

done



