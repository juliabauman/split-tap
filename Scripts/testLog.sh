#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --partition=larsms
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

#conda activate split-tap before running

ml load biology
ml load trim_galore/0.5.0
ml load samtools
ml load python/3.12.1

DATA_DIR="admera/test" #change if raw fastqs in differently named folder
WELLS=("A1")
FILE_FOLDERS=()

for well in "${WELLS[@]}"; do
    for folder in "${DATA_DIR}/${well}_"*; do  # Match folders like A1_1, A1_2, etc. - adjust folder structure as needed
        if [[ -d "$folder" ]]; then
            FILE_FOLDERS+=("$folder")
        fi
    done
done

#for index correction
I5_CDNA="TACTCCTT" #update as needed
I5_CROP="GGGGGGGG"

declare -A I7_INDEX_MAP
while IFS=',' read -r well sample_name i7_index; do
    well=$(echo "$well" | xargs)           # remove leading/trailing spaces
    i7_index=$(echo "$i7_index" | xargs)   # same for index
    I7_INDEX_MAP["$well"]="$i7_index"
done < i7_barcodes.csv


#mkdir "admera/test/barcoded_fastqs"

#Process all files thru barcoding
for folder in "${FILE_FOLDERS[@]}"; do
#
#    # Process cDNA files
#    if [[ "$folder" =~ _1$ ]]; then
#        files=($(find "$folder" -maxdepth 1 -type f -name "*.fastq.gz" | sort))
#        R1="${files[0]}"
#        R2="${files[1]}"
#        file_name=$(basename "$R1")
#        r2_file_name=$(basename "$R2")
#        sample=$(echo "$file_name" | cut -d'_' -f1)
#
#        if [[ -f "$R2" ]]; then
#            echo "Processing $R1 and $R2 from sample $sample"
#
#            mkdir -p "${folder}/trimmed"
#            mkdir -p "Data/barcoded_fastqs/${sample}"
#
#            # adapter trimming for cDNA  
#            trim_galore -o "${folder}/trimmed" --paired "${folder}/ind_corr/${file_name}" "${folder}/ind_corr/${r2_file_name}"
#            trimmed_R1="${folder}/trimmed/${file_name%.fastq.gz}_val_1.fq.gz"
#            trimmed_R2="${folder}/trimmed/${r2_file_name%.fastq.gz}_val_2.fq.gz"
#
#            #add in polyA trimming for R1, bc this fucks up STAR alignment if present
#            cutadapt -a "AAAAAAAAA" -o "${folder}/trimmed/R1_dbl_trim.fq.gz" "$trimmed_R1"
#
#            bash Scripts/parse_barcoding.sh "${folder}/trimmed/R1_dbl_trim.fq.gz" "$trimmed_R2"
#            mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/test_cDNA_barcode_head.fastq.gz"
#
#            #add index to the cell_barcode for each read
#            i7_index="${I7_INDEX_MAP[$sample]}"
#            python3 Scripts/combine_cell_bcs.py "Data/barcoded_fastqs/${sample}/test_cDNA_barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/test_cDNA_barcode_head_wIndex.fastq.gz" i7_index
#
#            echo "Finished processing $R1 and $R2"
#        else
#            echo "Warning: No matching R2 file found for $folder"
#        fi
#    fi
#
    # Process CROP folders

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
            #python Scripts/spacer_to_R1header.py "${R1}" "${R2}" "Data/fastq_files/${sample}/test_CROP_R1_spacerAdded.fq.gz" "Data/fastq_files/${sample}/test_CROP_R2_spacerAdded.fq.gz"
            
            #Then Parse barcoding on modified CROP folders & add index to CB
            #bash Scripts/parse_barcoding.sh "Data/fastq_files/${sample}/test_CROP_R1_spacerAdded.fq.gz" "Data/fastq_files/${sample}/test_CROP_R2_spacerAdded.fq.gz"
            #mv "Data/barcoded_fastqs/process/barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head.fastq.gz"
 
            i7_index="${I7_INDEX_MAP[$sample]}"
            echo "Annotating $sample reads with index $i7_index..."
            #python3 Scripts/combine_cell_bcs.py "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head.fastq.gz" "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_wIndex.fastq.gz" $i7_index

            #bash Scripts/spacer_to_header.sh "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_wIndex.fastq.gz" "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_spacer1.fastq.gz"
 
            #Then extract second spacer from barcoded R1 folder
            #python3 Scripts/extract_spacer2.py "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_spacer1.fastq.gz" "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_wSpacers.fastq.gz"
            echo "Finished processing $R1 and $R2 CROP files."
        else
            echo "Warning: No matching R2 file found for $R1"
        fi
    fi
done


#Combine all wells' files FOR cDNA & CROP, respectively
#zcat Data/barcoded_fastqs/*/cDNA_barcode_head_wIndex.fastq.gz | gzip > Data/barcoded_fastqs/all_cDNA_barcoded.fastq.gz
#zcat Data/barcoded_fastqs/*/CROP_barcode_head_wSpacers.fastq.gz | gzip > Data/barcoded_fastqs/all_CROP_barcoded.fastq.gz


#mkdir Data/BAM

##for cDNA folder:
#bash Scripts/star_align.sh Data/barcoded_fastqs/all_cDNA_barcoded.fastq.gz Data/BAM/
#samtools sort "Data/BAM/Aligned.out.bam" -o "Data/BAM/Aligned_sorted.out.bam"
#samtools index "Data/BAM/Aligned_sorted.out.bam"

ml load system
ml load ncurses/6.0
ml load py-biopython/1.79_py39
ml load py-pysam/0.18.0_py39

#Add necessary tags to the BAM folder (for UMI & CB)
#python3 Scripts/add_tags.py "Data/BAM/Aligned_sorted.out.bam" "Data/BAM/Aligned_sorted_wTags.out.bam"

#Feature counts & umi-tools grouping/dedup for cDNA
#bash Scripts/feat_cts.sh "Data/BAM/Aligned_sorted_wTags.out.bam" "Data/BAM" #now includes chimera correction

ml load py-pandas/2.2.1_py312

#Extract guides into txt folder, then assign guide pair id
#python3 Scripts/extract_guides.py "Data/barcoded_fastqs/${sample}/test_CROP_barcode_head_wSpacers.fastq.gz" "Data/CROP_summary.txt" "Data/CROP_summary_stats.txt" "Data/CROP_read_counts.txt"
#python3 Scripts/theseus_impl.py "Data/CROP_read_counts.txt" "Data/CROP_chimera_filt.txt" "0.25"

#filter CROP_summary by CROP_chimera_filt
echo "Loading CROP_filtered.txt into memory..."
awk '
    NR==FNR {
        key = $1 FS $2 FS $3 FS $4 FS $5;
        filter[key] = 1;
        next;
    }
    FNR==1 {print $0; next}
    { 
        key = $1 FS $2 FS $3 FS $4 FS $5;
        if (key in filter) print $0;  # If it exists in the hash table, print it
    }
' "Data/CROP_chimera_filt.txt" "Data/CROP_summary.txt" > "Data/CROP_summary_filt.txt"

echo "Filtering complete. Saved to CROP_summary_filt.txt"

python3 Scripts/id_guides.py "Data/CROP_summary_filt.txt" "all_spacer_pairs_rc.csv" "Data/CROP_ids_dedup.txt" "Data/CROP_counts.txt"

#Generate a file that calls the tRNA for the top guide for each cell (to use as covariate in sceptre; it will only analyze single-pert cells anyway)
python3 Scripts/trna_assign.py

#Combine cDNA & CROP count tables
#(cat Data/CROP_counts.tsv; sed 1d Data/cDNA_counts.txt) > Data/guide_cDNA_counts.tsv
#sed -i '1s/Guide1-Guide2_Name/gene/' Data/guide_cDNA_counts.tsv

