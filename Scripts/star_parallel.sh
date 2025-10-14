#!/bin/bash
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --job-name=star_align
#SBATCH --output=logs/star_%A_%a.out

ml load biology
ml load star/2.7.10b

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
  echo "❌ SLURM_ARRAY_TASK_ID is not set. Did you forget to use --array=?" >&2
  exit 1
fi

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cdna_input_files.txt | xargs)
sample_name=$(basename $(dirname "$sample"))
outdir="Data/star_outputs/$sample_name"

echo "➡️ SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
echo "➡️ sample = '$sample'"
echo "➡️ sample_name = '$sample_name'"

mkdir -p "$outdir"

#already generated genome

STAR --runThreadN 3 \
     --genomeDir STAR/full_transcripts \
     --readFilesIn "$sample" \
     --sjdbGTFfile STAR/full_transcripts.gtf \
     --outSAMmultNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "$outdir/" \
     --chimSegmentMin 20 \
     --limitBAMsortRAM 30000000000 \
     --outFilterMultimapNmax 30 \
     --chimOutType Junctions SeparateSAMold \
     --outFilterScoreMinOverLread .5 \
     --outFilterMatchNminOverLread .5 \
     --readFilesCommand zcat \
     --outReadsUnmapped Fastx
