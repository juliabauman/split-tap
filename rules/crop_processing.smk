# Rules for CROP/guide processing
# Different workflows for single vs paired guide modes

import math

# Helper to get CROP FASTQ files from folder
def get_crop_fastqs(wildcards):
    """Get R1 and R2 FASTQ files from CROP folder"""
    folder = WELL_FOLDERS["crop"][wildcards.well]
    fastqs = sorted(glob.glob(os.path.join(folder, "*.fastq.gz")) +
                    glob.glob(os.path.join(folder, "*.fq.gz")))
    if len(fastqs) < 2:
        raise ValueError(f"Expected at least 2 FASTQ files in {folder}, found {len(fastqs)}")
    return {"r1": fastqs[0], "r2": fastqs[1]}

# PAIRED MODE RULES
if GUIDE_MODE == "paired":

    rule crop_spacer_to_r1header: #split-pipe removes R2 file, so need to get guide data from there before extracting cell barcodes.
        """Extract spacer from R2 and add to R1 header"""
        input:
            unpack(get_crop_fastqs)
        output:
            r1 = temp(f"{OUTPUT_DIR}/fastq_files/{{well}}/CROP_R1_spacerAdded.fq.gz"),
            r2 = temp(f"{OUTPUT_DIR}/fastq_files/{{well}}/CROP_R2_spacerAdded.fq.gz")
        shell:
            """
            ml python/3.6
            python3 {SCRIPTS_DIR}/spacer_to_R1header.py {input.r1} {input.r2} {output.r1} {output.r2}
            """

    rule crop_parse_barcoding_paired:
        """Parse barcoding on spacer-modified CROP files"""
        input:
            r1 = f"{OUTPUT_DIR}/fastq_files/{{well}}/CROP_R1_spacerAdded.fq.gz",
            r2 = f"{OUTPUT_DIR}/fastq_files/{{well}}/CROP_R2_spacerAdded.fq.gz"
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head.fastq.gz")
        params:
            genome_dir = GENOME_DIR,
            output_dir = OUTPUT_DIR
        resources:
            mem_mb = lambda wildcards, attempt: int(config["resources"]["parse_mem"].replace("G", "000")) * attempt
        threads: config["threads"]["parse"]
        shell:
            """
            ml biology
            ml star
            bash {SCRIPTS_DIR}/parse_barcoding.sh {input.r1} {input.r2} {params.output_dir}/barcoded_fastqs
            mv {params.output_dir}/barcoded_fastqs/process/barcode_head.fastq.gz {output}
            """

    rule crop_add_index_paired:
        """Add i7 index to CROP cell barcodes"""
        input:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head.fastq.gz"
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wIndex.fastq.gz")
        params:
            i7_index = get_i7_index
        threads: config["threads"]["pigz"]
        shell:
            """
            bash {SCRIPTS_DIR}/combine_cell_bcs.sh {input} {output} {params.i7_index}
            """

    rule crop_spacer_to_header:
        """Extract first spacer from read sequence to header"""
        input:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wIndex.fastq.gz"
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_spacer1.fastq.gz")
        shell:
            """
            bash {SCRIPTS_DIR}/spacer_to_header.sh {input} {output}
            """

    rule crop_extract_spacer2:
        """Extract second spacer from read sequence to header"""
        input:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_spacer1.fastq.gz"
        output:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wSpacers.fastq.gz"
        shell:
            """
            ml python/3.6
            python3 {SCRIPTS_DIR}/extract_spacer2.py {input} {output}
            """

    rule crop_combine_wells_paired:
        """Combine all CROP files from all wells"""
        input:
            expand(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wSpacers.fastq.gz",
                   well=CROP_WELLS)
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/all_CROP_barcoded.fastq.gz")
        threads: config["threads"]["pigz"]
        shell:
            """
            zcat {input} | pigz -p {threads} > {output}
            """

    rule crop_extract_guides_paired:
        """Extract cell barcode, UMI, and spacer info to table"""
        input:
            f"{OUTPUT_DIR}/barcoded_fastqs/all_CROP_barcoded.fastq.gz"
        output:
            summary = f"{OUTPUT_DIR}/CROP_summary.txt",
            stats = f"{OUTPUT_DIR}/CROP_summary_stats.txt",
            read_counts = f"{OUTPUT_DIR}/CROP_read_counts.txt"
        shell:
            """
            ml python/3.6
            python3 {SCRIPTS_DIR}/extract_guides.py {input} {output.summary} {output.stats} {output.read_counts}
            """

    # Split into chunks for parallel UMI collapse
    checkpoint crop_split_for_collapse:
        """Split CROP read counts into chunks for parallel UMI collapse"""
        input:
            f"{OUTPUT_DIR}/CROP_read_counts.txt"
        output:
            directory(f"{OUTPUT_DIR}/CROP_reads_chunk")
        params:
            num_chunks = 58  # Same as array size in original
        shell:
            """
            mkdir -p {output}
            # Skip header and split into chunks
            tail -n +2 {input} > {output}/no_header.tmp
            split -n l/{params.num_chunks} -d -a 2 {output}/no_header.tmp {output}/chunk_
            # Add header back to each chunk and convert to TSV
            for chunk in {output}/chunk_*; do
                base=$(basename $chunk)
                (head -n 1 {input} && cat $chunk) > {output}/${{base}}.tsv
                rm $chunk
            done

            rm {output}/no_header.tmp
            """

    def get_chunk_files(wildcards):
        """Get all chunk files after checkpoint"""
        checkpoint_output = checkpoints.crop_split_for_collapse.get(**wildcards).output[0]
        chunk_files = glob.glob(os.path.join(checkpoint_output, "*.tsv"))
        return chunk_files

    rule crop_collapse_umis:
        """Collapse UMIs within each chunk"""
        input:
            f"{OUTPUT_DIR}/CROP_reads_chunk/chunk_{{chunk}}.tsv"
        output:
            collapsed = temp(f"{OUTPUT_DIR}/umi_collapse_temp/output_chunk_{{chunk}}.tsv"),
            umi_log = temp(f"{OUTPUT_DIR}/umi_collapse_temp/umi_corrections_chunk_{{chunk}}.tsv"),
            spacer_log = temp(f"{OUTPUT_DIR}/umi_collapse_temp/spacer_corrections_chunk_{{chunk}}.tsv")
        shell:
            """
            mkdir -p {OUTPUT_DIR}/umi_collapse_temp
            python3 {SCRIPTS_DIR}/collapse_umis.py {input} {output.collapsed} {output.umi_log} {output.spacer_log}
            """

    def get_all_collapsed_chunks(wildcards):
        """Get all collapsed chunk files"""
        checkpoint_output = checkpoints.crop_split_for_collapse.get(**wildcards).output[0]
        chunk_files = glob.glob(os.path.join(checkpoint_output, "*.tsv"))
        chunk_ids = [os.path.basename(f).replace("chunk_", "").replace(".tsv", "") for f in chunk_files]
        return expand(f"{OUTPUT_DIR}/umi_collapse_temp/output_chunk_{{chunk}}.tsv", chunk=chunk_ids)

    rule crop_combine_collapsed:
        """Combine all collapsed UMI chunks"""
        input:
            get_all_collapsed_chunks
        output:
            f"{OUTPUT_DIR}/umi_collapse_all.tsv"
        shell:
            """
            head -n 1 {input[0]} > {output}
            tail -n +2 -q {input} >> {output}
            """

    rule crop_theseus_paired:
        """Chimera correction using theseus"""
        input:
            f"{OUTPUT_DIR}/umi_collapse_all.tsv"
        output:
            f"{OUTPUT_DIR}/CROP_chimera_corr_filt.txt"
        params:
            threshold = config["chimera_threshold"]
        shell:
            """
            python3 {SCRIPTS_DIR}/theseus_impl_paired.py {input} {output} {params.threshold} {OUTPUT_DIR}/tpt_histogram.pdf
            """

    rule crop_filter_by_theseus:
        """Filter CROP summary by theseus results"""
        input:
            filter_file = f"{OUTPUT_DIR}/CROP_chimera_corr_filt.txt",
            summary = f"{OUTPUT_DIR}/CROP_summary.txt"
        output:
            f"{OUTPUT_DIR}/CROP_summary_corr_filt.txt"
        shell:
            """
            awk -F'\t' 'BEGIN{{OFS="\t"}}
            # First file: build filter set from CROP_chimera_corr_filt
            NR==FNR {{
                if (FNR==1) {{
                # Identify columns by name
                    for (i=1;i<=NF;i++) {{
                        if ($i=="spacer1_spacer2") s_col=i;
                        else if ($i=="cell") c_col=i;
                        else if ($i=="umi") u_col=i;
                    }}
                    next;
                }}
                key = $s_col OFS $c_col OFS $u_col;
                filter[key]=1;
                next;
            }}
            # Second file: CROP_summary.txt
            FNR==1 {{
                # Identify columns in summary
                for (i=1;i<=NF;i++) {{
                    if ($i=="spacer1_spacer2") s2=i;
                    else if ($i=="cell") c2=i;
                    else if ($i=="umi") u2=i;
                }}
                print;  # print header as-is
                next;
            }}
            {{
                key = $s2 OFS $c2 OFS $u2;
                if (key in filter) print;
            }}          
            ' {input.filter_file} {input.summary} > {output}
            """

    rule crop_id_guides_paired:
        """Assign guide pair IDs and generate counts"""
        input:
            summary = f"{OUTPUT_DIR}/CROP_summary_corr_filt.txt",
            spacer_pairs = config["spacer_pairs_csv"]
        output:
            ids = f"{OUTPUT_DIR}/CROP_ids_dedup.txt",
            counts = f"{OUTPUT_DIR}/CROP_counts.txt"
        shell:
            """
            python3 {SCRIPTS_DIR}/id_guides.py {input.summary} {input.spacer_pairs} {output.ids} {output.counts}
            """

    rule crop_trna_assign:
        """Assign tRNA for top guide in each cell"""
        input:
            f"{OUTPUT_DIR}/CROP_ids_dedup.txt"
        output:
            f"{OUTPUT_DIR}/tRNA_assignments.txt"
        shell:
            """
            python3 {SCRIPTS_DIR}/trna_assign.py {input} {output}
            """

# SINGLE MODE RULES
else:  # GUIDE_MODE == "single"

    rule crop_parse_barcoding_single:
        """Parse barcoding on CROP files"""
        input:
            unpack(get_crop_fastqs)
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head.fastq.gz")
        params:
            genome_dir = GENOME_DIR,
            output_dir = OUTPUT_DIR
        resources:
            mem_mb = lambda wildcards, attempt: int(config["resources"]["parse_mem"].replace("G", "000")) * attempt
        threads: config["threads"]["parse"]
        shell:
            """
            ml biology
            ml star
            bash {SCRIPTS_DIR}/parse_barcoding.sh {input.r1} {input.r2} {params.output_dir}/barcoded_fastqs
            mv {params.output_dir}/barcoded_fastqs/process/barcode_head.fastq.gz {output}
            """

    rule crop_add_index_single:
        """Add i7 index to CROP cell barcodes"""
        input:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head.fastq.gz"
        output:
            f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wIndex.fastq.gz"
        params:
            i7_index = get_i7_index
        threads: config["threads"]["pigz"]
        shell:
            """
            bash {SCRIPTS_DIR}/combine_cell_bcs.sh {input} {output} {params.i7_index}
            """

    rule crop_combine_wells_single:
        """Combine all CROP files from all wells"""
        input:
            expand(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/CROP_barcode_head_wIndex.fastq.gz",
                   well=CROP_WELLS)
        output:
            temp(f"{OUTPUT_DIR}/barcoded_fastqs/all_CROP_barcoded.fastq.gz")
        threads: config["threads"]["pigz"]
        shell:
            """
            zcat {input} | pigz -p {threads} > {output}
            """

    rule crop_process_guides_single:
        """Extract single guides from CROP file"""
        input:
            fastq = f"{OUTPUT_DIR}/barcoded_fastqs/all_CROP_barcoded.fastq.gz",
            guide_lib = config["guide_library_csv"]
        output:
            f"{OUTPUT_DIR}/guides.txt"
        shell:
            """
            ml biology
            ml py-biopython
            python3 {SCRIPTS_DIR}/process_guides_new.py {input.fastq} {output} {input.guide_lib}
            """

    rule crop_theseus_single:
        """Chimera correction using theseus"""
        input:
            f"{OUTPUT_DIR}/guides.txt"
        output:
            f"{OUTPUT_DIR}/guides_filtered.txt"
        params:
            threshold = config["chimera_threshold"]
        resources:
            mem_mb=95000
        shell:
            """
            python3 {SCRIPTS_DIR}/theseus_impl.py {input} {output} {params.threshold} {OUTPUT_DIR}/tpt_histogram.pdf
            """

    rule crop_count_guides_single:
        """Count guides after filtering"""
        input:
            f"{OUTPUT_DIR}/guides_filtered.txt"
        output:
            f"{OUTPUT_DIR}/guide_counts.txt"
        shell:
            """
            python3 {SCRIPTS_DIR}/count_guides.py {input} {output}
            """
