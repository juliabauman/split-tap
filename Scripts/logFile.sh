#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=100G
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8

#conda activate split-tap before running

ml load biology
ml load samtools
ml load htslib


DATA_DIR="admera" #change if raw fastqs in differently named folder
#WELLS=("A1" "A3" "A4" "A5" "A6" "A7" "A8" "A9" "A10" "A11" "A12" "B1" "B2" "B3" "B4" "B5" "B6" "B7" "B8" "B9" "B10" "B11" "B12" "C2" "C3" "C4")
WELLS=("C5" "C6" "C7" "C8") 
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


#mkdir -p Data/barcoded_fastqs
#mkdir -p Data/fastq_files


#Process all files thru barcoding
for folder in "${FILE_FOLDERS[@]}"; do

    # Process cDNA files
    if [[ "$folder" =~ _1$ ]]; then
        files=($(find "$folder" -maxdepth 1 -type f -name "*.fastq.gz" | sort))
        R1="${files[0]}"
        R2="${files[1]}"
        file_name=$(basename "$R1")
        r2_file_name=$(basename "$R2")
        sample=$(echo "$file_name" | cut -d'_' -f1)

        if [[ -f "$R2" ]]; then
            echo "Processing $R1 and $R2 from sample $sample"

            correct_i7="${I7_INDEX_MAP[$sample]}"

            #mkdir -p "${folder}/trimmed"
            mkdir -p "Data/barcoded_fastqs/${sample}"

            bash Scripts/parse_barcoding.sh "${R1}" "${R2}"
            mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/cDNA_barcode_head.fastq.gz"

            # adapter trimming for cDNA, add in polyA trimming for R1, bc this fucks up STAR alignment if present
            cutadapt -j 20 \
              -a "AGATCGGAAGAGC" \
              -a A{10}N{100} \
              -q 0,20 \
              --minimum-length 20 \
              -o "Data/barcoded_fastqs/${sample}/cDNA_barcode_head_trimmed.fastq.gz" \
              "Data/barcoded_fastqs/${sample}/cDNA_barcode_head.fastq.gz"

            #add index to the cell_barcode for each read
            bash Scripts/combine_cell_bcs.sh "Data/barcoded_fastqs/${sample}/cDNA_barcode_head_trimmed.fastq.gz" "Data/barcoded_fastqs/${sample}/cDNA_barcode_head_wIndex.fastq.gz" ${correct_i7}

            echo "Finished processing $R1 and $R2"
        else
            echo "Warning: No matching R2 file found for $folder"
        fi
    fi


   #Process CROP files

    if [[ "$folder" =~ _2$ ]]; then
        files=($(find "$folder" -maxdepth 1 -type f -name "*.fastq.gz" | sort))
        R1="${files[0]}"
        R2="${files[1]}"
        file_name=$(basename "$R1")
        r2_file_name=$(basename "$R2")
        sample=$(echo "$file_name" | cut -d'_' -f1)

        # Ensure the second file exists before processing
        if [[ -f "$R2" ]]; then
            echo "Processing $R1 and $R2 from sample $sample"

            #mkdir -p "Data/fastq_files/${sample}"

            #Move CROP R2-contained spacer seqs to R1 file read header (bc Parse pipeline removes R2)
            python Scripts/spacer_to_R1header.py "${R1}" "${R2}" "Data/fastq_files/${sample}/test_CROP_R1_spacerAdded.fq.gz" "Data/fastq_files/${sample}/test_CROP_R2_spacerAdded.fq.gz"

            #Then Parse barcoding on modified CROP folders & add index to CB
            bash Scripts/parse_barcoding.sh "Data/fastq_files/${sample}/test_CROP_R1_spacerAdded.fq.gz" "Data/fastq_files/${sample}/test_CROP_R2_spacerAdded.fq.gz"
            mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head.fastq.gz"

            i7_index="${I7_INDEX_MAP[$sample]}"
            echo "Annotating $sample reads with index $i7_index..."
            bash Scripts/combine_cell_bcs.sh "Data/barcoded_fastqs/${sample}/CROP_barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head_wIndex.fastq.gz" $i7_index

            bash Scripts/spacer_to_header.sh "Data/barcoded_fastqs/${sample}/CROP_barcode_head_wIndex.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head_spacer1.fastq.gz"

            #Then extract second spacer from barcoded R1 folder
            python3 Scripts/extract_spacer2.py "Data/barcoded_fastqs/${sample}/CROP_barcode_head_spacer1.fastq.gz" "Data/barcoded_fastqs/${sample}/CROP_barcode_head_wSpacers.fastq.gz"
            echo "Finished processing $R1 and $R2 CROP files."
        else
            echo "Warning: No matching R2 file found for $R1"
        fi
    fi
done

ml load system
ml load pigz/2.4
ml load py-pysam
ml load py-pandas

#Combine all wells' files for CROP (not cDNA - too big, need to parallelize)
#zcat Data/barcoded_fastqs/*/CROP_barcode_head_wSpacers.fastq.gz | pigz -p 20 > Data/barcoded_fastqs/all_CROP_barcoded.fastq.gz

#Extract guides into txt folder, then assign guide pair id
#python3 Scripts/extract_guides.py "Data/barcoded_fastqs/all_CROP_barcoded.fastq.gz" "Data/CROP_summary.txt" "Data/CROP_summary_stats.txt" "Data/CROP_read_counts.txt"
#sbatch Scripts/collapse_batch.sh ## Need to work in a way to prevent future steps from running before this one is done
#head -n 1 Data/umi_collapse_temp/output_chunk_aa.tsv > Data/umi_collapse_all.tsv
#for f in Data/umi_collapse_temp/output_*.tsv; do
#    tail -n +2 "$f" >> Data/umi_collapse_all.tsv
#done
#python3 Scripts/theseus_impl.py "Data/umi_collapse_all.tsv" "Data/CROP_chimera_corr_filt.txt" "0.25"


##Parallelize STAR alignment, tag addition & featureCounts for cDNA
#ls Data/barcoded_fastqs/*/cDNA_barcode_head_wIndex.fastq.gz > cdna_input_files.txt
#jobid=$(sbatch --wait --array=1-26 Scripts/star_parallel.sh)
#echo "Submitted STAR array job with ID $jobid"

#ls Data/star_outputs/*/Aligned.sortedByCoord.out.bam > cdna_input_bams.txt
#tag_jobid=$(sbatch --wait Scripts/tag_parallel.sh)

#ls Data/star_outputs/*/Aligned.sortedByCoord.out_tagged.bam > aligned_tagged_bams.txt
#mkdir -p Data/final_outputs
ft_cts_jobid=$(sbatch --wait Scripts/cdna_summ.sh | awk '{print $4}')
echo "Submitted tagging job with ID $ft_cts_jobid"

##Combine all samples' cDNA_counts files
head -n 1 Data/final_outputs/*_cDNA_counts.tsv | head -n 1 > Data/final_outputs/combined_cDNA_counts.tsv
tail -n +2 -q Data/final_outputs/*_cDNA_counts.tsv >> Data/combined_cDNA_counts.tsv


##Finish processing CROP
awk '
    NR==FNR {
        key = $1 FS $2 FS $4 FS $5;
        filter[key] = 1;
        next;
    }
    FNR==1 { print; next }  # print header
    {
        key = $3 FS $4 FS $1 FS $2;
        if (key in filter) print;
    }
' Data/CROP_chimera_corr_filt.txt Data/CROP_summary.txt > Data/CROP_summary_corr_filt.txt

echo "Filtering complete. Saved to CROP_summary_filt.txt"

python3 Scripts/id_guides.py "Data/CROP_summary_corr_filt.txt" "all_spacer_pairs_rc.csv" "Data/CROP_ids_dedup.txt" "Data/final_outputs/CROP_counts.txt"

#Generate a file that calls the tRNA for the top guide for each cell (to use as covariate in sceptre; it will only analyze single-pert cells anyway)
python3 Scripts/trna_assign.py


#Combine cDNA & CROP count tables - this WON'T work as written, columns need to be reordered!!
#(cat Data/CROP_counts.tsv; sed 1d Data/cDNA_counts.txt) > Data/guide_cDNA_counts.tsv
#sed -i '1s/Guide1-Guide2_Name/gene/' Data/guide_cDNA_counts.tsv

