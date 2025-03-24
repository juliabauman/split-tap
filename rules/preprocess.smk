rule split_pipe:
    input:
        r1=config["fastq"]["r1"],
        r2=config["fastq"]["r2"]
    output:
        directory(f"{config['output_dir']}/barcoded_fastqs")
    shell:
        """
        split-pipe --mode pre --chemistry v1 \
            --genome_dir {config[genome_dir]} \
            --fq1 {input.r1} --fq2 {input.r2} --output_dir {output}
        """

rule star_align:
    input:
        f"{config['output_dir']}/barcoded_fastqs/process/barcode_head.fastq.gz"
    output:
        bam=f"{config['output_dir']}/BAM/Aligned.out.bam",
        sorted_bam=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam",
        bai=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam.bai"
    shell:
        """
        ml load biology
        ml load star/2.7.10b
        ml load samtools/1.16.1

        STAR --runThreadN {config[threads]} \
            --genomeDir {config[star_index]} \
            --readFilesIn {input} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {config[output_dir]}/BAM/

        samtools sort {output.bam} -o {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

