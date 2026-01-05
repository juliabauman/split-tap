# Rules for STAR alignment, tagging, and featureCounts
# Workflow: STAR align -> Add tags (CB, UB) -> featureCounts -> umi_tools count

rule star_align:
    """Align cDNA reads with STAR"""
    input:
        f"{OUTPUT_DIR}/barcoded_fastqs/{{well}}/cDNA_barcode_head_wIndex.fastq.gz"
    output:
        bam = f"{OUTPUT_DIR}/star_outputs/{{well}}/Aligned.sortedByCoord.out.bam",
        log_final = f"{OUTPUT_DIR}/star_outputs/{{well}}/Log.final.out",
        log_out = f"{OUTPUT_DIR}/star_outputs/{{well}}/Log.out",
        log_progress = f"{OUTPUT_DIR}/star_outputs/{{well}}/Log.progress.out"
    params:
        genome_dir = STAR_GENOME_DIR,
        gtf = STAR_GTF,
        outdir = lambda wildcards: f"{OUTPUT_DIR}/star_outputs/{wildcards.well}"
    threads: config["threads"]["star"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["resources"]["star_mem"].replace("G", "000")) * attempt
    shell:
        """
        ml biology
        ml star/2.7.10b
        mkdir -p {params.outdir}

        STAR --runThreadN {threads} \
             --genomeDir {params.genome_dir} \
             --readFilesIn {input} \
             --sjdbGTFfile {params.gtf} \
             --outSAMmultNmax 1 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.outdir}/ \
             --chimSegmentMin 20 \
             --limitBAMsortRAM 30000000000 \
             --outFilterMultimapNmax 30 \
             --chimOutType Junctions SeparateSAMold \
             --outFilterScoreMinOverLread .5 \
             --outFilterMatchNminOverLread .5 \
             --readFilesCommand zcat \
             --outReadsUnmapped Fastx
        """

rule add_tags:
    """Add CB and UB tags to BAM file"""
    input:
        f"{OUTPUT_DIR}/star_outputs/{{well}}/Aligned.sortedByCoord.out.bam"
    output:
        bam = f"{OUTPUT_DIR}/star_outputs/{{well}}/Aligned.sortedByCoord.out_tagged.bam",
        bai = f"{OUTPUT_DIR}/star_outputs/{{well}}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ml biology
        ml samtools
        ml python/3.6
        ml py-pysam
        samtools index {input}
        python3 {SCRIPTS_DIR}/add_tags.py {input} {output.bam}
        """

rule featurecounts_and_count:
    """Run featureCounts and count UMIs with umi_tools"""
    input:
        bam = f"{OUTPUT_DIR}/star_outputs/{{well}}/Aligned.sortedByCoord.out_tagged.bam"
    output:
        counts = f"{OUTPUT_DIR}/final_outputs/{{well}}_cDNA_counts.tsv",
        feature_counts = f"{OUTPUT_DIR}/featureCounts/{{well}}/amplicon_assigned_{{well}}.txt",
        feature_summary = f"{OUTPUT_DIR}/featureCounts/{{well}}/amplicon_assigned_{{well}}.txt.summary",
        tagged_bam = f"{OUTPUT_DIR}/featureCounts/{{well}}/Aligned.sortedByCoord.out_tagged.bam.featureCounts.bam",
        filtered_bam = temp(f"{OUTPUT_DIR}/featureCounts/{{well}}/Filtered_{{well}}_tagged.bam"),
        filtered_sorted = f"{OUTPUT_DIR}/featureCounts/{{well}}/Filtered_{{well}}_tagged_sorted.bam",
        filtered_bai = f"{OUTPUT_DIR}/featureCounts/{{well}}/Filtered_{{well}}_tagged_sorted.bam.bai"
    params:
        gtf = STAR_GTF,
        outdir = lambda wildcards: f"{OUTPUT_DIR}/featureCounts/{wildcards.well}"
    threads: config["threads"]["featurecounts"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["resources"]["featurecounts_mem"].replace("G", "000")) * attempt
    shell:
        """
        mkdir -p {params.outdir}
        ml biology
        ml samtools

        # Run featureCounts
        featureCounts -a {params.gtf} \
            -o {output.feature_counts} \
            -g gene_id -R BAM -M --primary -T {threads} {input.bam}

        # Summarize UMI–CB–gene
        samtools view {output.tagged_bam} | \
        awk '
          {{
            umi=""; cell=""; gene="";
            for(i=1;i<=NF;i++) {{
              if ($i ~ /^UB:Z:/) umi = substr($i, 6);
              if ($i ~ /^CB:Z:/) cell = substr($i, 6);
              if ($i ~ /^XT:Z:/) gene = substr($i, 6);
            }}
            if (umi && cell && gene)
              print umi "\\t" cell "\\t" gene;
          }}
        ' | sort | uniq -c | sort -nr | \
        awk 'BEGIN {{ print "Num_Obs\\tMolecular_Barcode\\tCell_Barcode\\tGene" }}
             NF==4 {{ printf "%s\\t%s\\t%s\\t%s\\n", $1, $2, $3, $4 }}' > {params.outdir}/umi_gene_cell_raw_counts.txt

        # Sort and index the filtered BAM (placeholder - using tagged BAM for now)
        cp {output.tagged_bam} {output.filtered_bam}
        samtools sort {output.filtered_bam} -o {output.filtered_sorted}
        samtools index {output.filtered_sorted}

        # Final UMI counting with umi_tools
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
                        --extract-umi-method=tag --umi-tag=UB \
                        --per-cell --cell-tag=CB \
                        -I {output.filtered_sorted} \
                        -S {output.counts}
        """
