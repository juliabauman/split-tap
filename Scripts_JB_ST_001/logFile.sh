#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8

#conda activate split-tap before running

ml load biology
ml load samtools
ml load htslib


DATA_DIR="../screen/usftp21.novogene.com/01.RawData" #change if raw fastqs in differently named folder
WELLS=("a1" "a2" "a3" "a4" "a5" "a6" "a7" "b1" "b2" "b3" "b4" "b5" "b7" "c1" "c2" "c3" "c4" "c5" "c6" "c7")
FILE_FOLDERS=()

for well in "${WELLS[@]}"; do
    for folder in "${DATA_DIR}/${well}_"*; do  # Match folders like A1_1, A1_2, etc. - adjust folder structure as needed
        if [[ -d "$folder" ]]; then
            FILE_FOLDERS+=("$folder")
        fi
    done
done


declare -A I7_INDEX_MAP
while IFS=',' read -r well sample_name i7_index; do
    well=$(echo "$well" | xargs)           # remove leading/trailing spaces
    i7_index=$(echo "$i7_index" | xargs)   # same for index
    I7_INDEX_MAP["$well"]="$i7_index"
done < i7_barcodes.csv


mkdir -p Data/barcoded_fastqs
mkdir -p Data/fastq_files


#Process all files thru barcoding
#for folder in "${FILE_FOLDERS[@]}"; do

    # Process cDNA files
   # if [[ "$folder" =~ _1$ ]]; then
   #     files=($(find "$folder" -maxdepth 1 -type f -name "*.fq.gz" | sort))
   #     R1="${files[0]}"
   #     R2="${files[1]}"
   #     file_name=$(basename "$R1")
   #     r2_file_name=$(basename "$R2")
   #     sample=$(echo "$file_name" | cut -d'_' -f1)

   #     if [[ -f "$R2" ]]; then
   #         echo "Processing $R1 and $R2 from sample $sample"

   #         correct_i7="${I7_INDEX_MAP[$sample]}"

   #         mkdir -p "${folder}/trimmed"
   #         mkdir -p "Data/barcoded_fastqs/${sample}"

   #         # adapter trimming for cDNA, add in polyA trimming for R1, bc this fucks up STAR alignment if present
   #         cutadapt -j 20 \
   #           -a "AGATCGGAAGAGC" -A "AGATCGGAAGAGC" \
   #           -a "AAAAAAAAAA" -A "AAAAAAAAAA" \
   #           -q 0,20 \
   #           --minimum-length 20 \
   #           -o "${folder}/trimmed/R1_trimmed.fq.gz" \
   #           -p "${folder}/trimmed/R2_trimmed.fq.gz" \
   #           "$R1" "$R2"

   #         bash Scripts/parse_barcoding.sh "${folder}/trimmed/R1_trimmed.fq.gz" "${folder}/trimmed/R2_trimmed.fq.gz"
   #         mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/cDNA_barcode_head.fastq.gz"

   #         #add index to the cell_barcode for each read
   #         bash Scripts/combine_cell_bcs.sh "Data/barcoded_fastqs/${sample}/cDNA_barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/cDNA_barcode_head_wIndex.fastq.gz" ${correct_i7}

   #         echo "Finished processing $R1 and $R2"
   #     else
   #         echo "Warning: No matching R2 file found for $folder"
   #     fi
   # fi


   #Process CROP files

   # if [[ "$folder" =~ _2$ ]]; then
   #     files=($(find "$folder" -maxdepth 1 -type f -name "*.fq.gz" | sort))
   #     R1="${files[0]}"
   #     R2="${files[1]}"
   #     file_name=$(basename "$R1")
   #     r2_file_name=$(basename "$R2")
   #     sample=$(echo "$file_name" | cut -d'_' -f1)

   #     # Ensure the second file exists before processing
   #     if [[ -f "$R2" ]]; then
   #         echo "Processing $R1 and $R2 from sample $sample"

   #         #Then Parse barcoding on modified CROP folders & add index to CB
   #         bash Scripts/parse_barcoding.sh "${R1}" "${R2}"
   #         mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head.fastq.gz"

   #         i7_index="${I7_INDEX_MAP[$sample]}"
   #         echo "Annotating $sample reads with index $i7_index..."
   #         bash Scripts/combine_cell_bcs.sh "Data/barcoded_fastqs/${sample}/CROP_barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head_wIndex.fastq.gz" $i7_index

   #         echo "Finished processing $R1 and $R2 CROP files."
   #     else
   #         echo "Warning: No matching R2 file found for $R1"
   #     fi
   # fi
#done

ml load system
ml load pigz/2.4
ml load py-pysam
ml load py-pandas

#Combine all wells' files for CROP (not cDNA - too big, need to parallelize future steps)
#zcat Data/barcoded_fastqs/*/CROP_barcode_head_wIndex.fastq.gz | pigz -p 20 > Data/barcoded_fastqs/all_CROP_barcoded.fastq.gz

#Spacer extraction from CROP file
bash Scripts/runProcess_guides.sh

##Parallelize STAR alignment, tag addition & featureCounts for cDNA
#ls Data/barcoded_fastqs/*/cDNA_barcode_head_wIndex.fastq.gz > cdna_input_files.txt
#jobid=$(sbatch --array=1-26 Scripts/star_parallel.sh | awk '{print $4}')
#echo "Submitted STAR array job with ID $jobid"

#ls Data/star_outputs/*/Aligned.sortedByCoord.out.bam > cdna_input_bams.txt
#tag_jobid=$(sbatch --dependency=afterok:$jobid Scripts/tag_parallel.sh | awk '{print $4}')
#echo "Submitted tagging job with ID $tag_jobid"

#ls Data/star_outputs/*/Aligned.sortedByCoord_tagged.bam > aligned_tagged_bams.txt
#mkdir -p Data/final_outputs
#ft_cts_jobid=$(sbatch --dependency=afterok:$tag_jobid Scripts/cdna_summ.sh | awk '{print $4}')
#echo "Submitted tagging job with ID $ft_cts_jobid"

#Combine all samples' cDNA_counts files
#head -n 1 Data/final_outputs/*_cDNA_counts.tsv | head -n 1 > Data/final_outputs/combined_cDNA_counts.tsv
#tail -n +2 -q Data/final_outputs/*_cDNA_counts.tsv >> Data/final_outputs/combined_cDNA_counts.tsv


##Finish processing CROP


#Combine cDNA & CROP count tables
(cat Data/CROP_counts.tsv; sed 1d Data/cDNA_counts.txt) > Data/guide_cDNA_counts.tsv
sed -i '1s/Guide1-Guide2_Name/gene/' Data/guide_cDNA_counts.tsv

