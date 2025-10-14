#!/bin/bash
#SBATCH --job-name=ft_cts_thru_umi_tools
#SBATCH --output=logs/sample_%A_%a.log
#SBATCH --error=logs/sample_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-19

# Load modules
#ml load biology
#ml load samtools/1.16.1
#ml load ncurses/6.0
#ml load py-biopython/1.79_py39
#ml load py-pandas/2.2.1_py312

# Read input BAM path
BAM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" aligned_tagged_bams.txt)
SAMPLE=$(basename $(dirname "$BAM"))
BASENAME=$(basename "$BAM")

OUTDIR="Data/featureCounts/$SAMPLE"
mkdir -p "$OUTDIR"

# 1. Run featureCounts and tag reads
#featureCounts -a STAR/full_transcripts.gtf \
#  -o "$OUTDIR/amplicon_assigned_${SAMPLE}.txt" \
#  -g gene_name -R BAM -M --primary -T 4 "$BAM"

TAGGED_BAM="${OUTDIR}/${BASENAME}.featureCounts.bam"

# 2. Summarize UMI–CB–gene
#samtools view "$TAGGED_BAM" | \
#awk '
#  {
#    umi=""; cell=""; gene="";
#    for(i=1;i<=NF;i++) {
#      if ($i ~ /^UB:Z:/) umi = substr($i, 6);
#      if ($i ~ /^CB:Z:/) cell = substr($i, 6);
#      if ($i ~ /^XT:Z:/) gene = substr($i, 6);
#    }
#    if (umi && cell && gene)
#      print umi "\t" cell "\t" gene;
#  }
#' | sort | uniq -c | sort -nr | \
#awk 'BEGIN { print "Num_Obs\tMolecular_Barcode\tCell_Barcode\tGene" }
#     NF==4 { printf "%s\t%s\t%s\t%s\n", $1, $2, $3, $4 }' > "$OUTDIR/umi_gene_cell_raw_counts.txt"

# 3. Chimera filtering
#python3 Scripts/theseus_impl.py \
#  "$OUTDIR/umi_gene_cell_raw_counts.txt" \
#  "$OUTDIR/chimera_corrected.txt" 0.7
#
## 4. Filter BAM to keep "true" reads
#awk 'NR > 1 {print $2, $3, $4}' "$OUTDIR/chimera_corrected.txt" \
#  > "$OUTDIR/valid_keys.txt"

#samtools view -h "${TAGGED_BAM}" | \
#awk -v keyfile="$OUTDIR/valid_keys.txt" '
#  BEGIN {
#    OFS="\t";
#    while ((getline < keyfile) > 0) key[$1, $2, $3] = 1;
#  }
#  /^@/ { print; next }
#  {
#    umi=""; cell=""; gene="";
#    for(i=1;i<=NF;i++) {
#      if ($i ~ /^UB:Z:/) umi=substr($i,6);
#      if ($i ~ /^CB:Z:/) cell=substr($i,6);
#      if ($i ~ /^XT:Z:/) gene=substr($i,6);
#    }
#    if ((umi, cell, gene) in key) print $0;
#  }' | samtools view -b -o "$OUTDIR/Filtered_${SAMPLE}_tagged.bam" -
#
#samtools sort "$OUTDIR/Filtered_${SAMPLE}_tagged.bam" -o "$OUTDIR/Filtered_${SAMPLE}_tagged_sorted.bam"
#samtools index "$OUTDIR/Filtered_${SAMPLE}_tagged_sorted.bam"
#
#have to remove all modules to get umi_tools to work - maybe still won't work
#module purge 

# 5. Final UMI counting
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
                --extract-umi-method=tag --umi-tag=UB \
                --per-cell --cell-tag=CB \
                -I "$OUTDIR/Filtered_${SAMPLE}_tagged_sorted.bam" \
                -S Data/final_outputs/${SAMPLE}_cDNA_counts.tsv

