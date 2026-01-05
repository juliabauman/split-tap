# Rules for cDNA processing
# Workflow: Parse barcoding -> Cutadapt trimming -> Add i7 index

# Helper to get FASTQ files from folder
def get_cdna_fastqs(wildcards):
    """Get R1 and R2 FASTQ files from cDNA folder"""
    folder = WELL_FOLDERS["cdna"][wildcards.well]
    fastqs = sorted(glob.glob(os.path.join(folder, "*.fastq.gz")) +
                    glob.glob(os.path.join(folder, "*.fq.gz")))
    if len(fastqs) < 2:
        raise ValueError(f"Expected at least 2 FASTQ files in {folder}, found {len(fastqs)}")
    return {"r1": fastqs[0], "r2": fastqs[1]}

# Helper to get i7 index for a well
def get_i7_index(wildcards):
    """Read i7 index from CSV file"""
    import pandas as pd
    df = pd.read_csv(config["i7_barcodes_csv"])
    df.columns = df.columns.str.strip()
    df['well'] = df['well'].str.strip()
    well_data = df[df['well'] == wildcards.well]
    if len(well_data) == 0:
        raise ValueError(f"No i7 index found for well {wildcards.well}")
    return well_data['i7_index'].values[0].strip()

rule cdna_parse_barcoding:
    input:
        unpack(get_cdna_fastqs)
    output:
        temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head.fastq.gz")
    params:
        genome_dir = GENOME_DIR,
        output_dir = OUTPUT_DIR
    resources:
        mem_mb = lambda wildcards, attempt: int(config["resources"]["parse_mem"].replace("G", "000")) * attempt
    threads: config["threads"]["parse"]
    shell:
        """
        ml biology
        ml samtools
        split-pipe \
            --mode pre \
            --chemistry v1 \
            --genome_dir {params.genome_dir} \
            --fq1 {input.r1} \
            --fq2 {input.r2} \
            --output_dir {params.output_dir}/barcoded_fastqs

        mv {params.output_dir}/barcoded_fastqs/process/barcode_head.fastq.gz {output}
        """

rule cdna_cutadapt_trim:
    input:
        f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head.fastq.gz"
    output:
        temp(f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head_trimmed.fastq.gz")
    threads: config["threads"]["cutadapt"]
    shell:
        """
        ml biology
        ml py-cutadapt/1.18_py36
        ml python/3.6
        cutadapt -j {threads} \
            -a "AGATCGGAAGAGC" \
            -a A{{10}}N{{100}} \
            -q 0,20 \
            --minimum-length 20 \
            -o {output} \
            {input}
        """

rule cdna_add_index:
    input:
        f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head_trimmed.fastq.gz"
    output:
        f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head_wIndex.fastq.gz"
    params:
        i7_index = get_i7_index
    threads: config["threads"]["pigz"]
    shell:
        """
        bash {SCRIPTS_DIR}/combine_cell_bcs.sh {input} {output} {params.i7_index}
        """
