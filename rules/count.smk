rule feature_counts:
    input:
        bam=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam"
    output:
        txt=f"{config['output_dir']}/summary.txt",
        annotated=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam.featureCounts.bam",
        bai=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam.featureCounts.bam.bai"
    shell:
        """
        ml load biology
        ml load samtools/1.16.1
        featureCounts -t exon -g gene_id \
            -a {config[gtf]} \
            -o {output.txt} -R BAM {input.bam}
        mv "../Data/Aligned_sorted.out.bam.featureCounts.bam" {output.annotated}
        samtools sort {output.annotated} -o {output.annotated}.tmp
        mv {output.annotated}.tmp {output.annotated}
        samtools index {output.annotated}
        """

rule add_tags:
    input:
        bam=f"{config['output_dir']}/BAM/Aligned_sorted.out.bam.featureCounts.bam"
    output:
        tagged=f"{config['output_dir']}/BAM/Aligned_sorted_ftct_wTags.bam",
        bai=f"{config['output_dir']}/BAM/Aligned_sorted_ftct_wTags.bam.bai"
    shell:
        """
        ml load biology
        ml load samtools
        python ../Scripts/add_tags.py {input.bam} {output.tagged}
        samtools index {output.tagged}
        """

rule umi_tools_count:
    input:
        f"{config['output_dir']}/BAM/Aligned_sorted_ftct_wTags.bam"
    output:
        f"{config['output_dir']}/cdna_counts.tsv.gz"
    shell:
        """
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
            --extract-umi-method=tag --umi-tag=UB \
            -I {input} -S {output}
        """
