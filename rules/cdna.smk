from pathlib import Path


def _cdna_prepare_dir(sample):
    return OUTPUT_DIR / "intermediate" / sample / "cdna"


def cdna_prepared_fastq(sample, read):
    return _cdna_prepare_dir(sample) / f"{sample}_cdna_R{read}_prepared.fastq.gz"


def cdna_barcoded_fastq(sample):
    return OUTPUT_DIR / "barcoded_fastqs" / sample / "cDNA_barcode_head.fastq.gz"


def cdna_indexed_fastq(sample):
    return OUTPUT_DIR / "barcoded_fastqs" / sample / "cDNA_barcode_head_wIndex.fastq.gz"


def cdna_star_dir(sample):
    return OUTPUT_DIR / "cdna" / sample


def cdna_bam_paths(sample):
    base_dir = cdna_star_dir(sample)
    return {
        "bam": base_dir / "Aligned.out.bam",
        "sorted": base_dir / "Aligned_sorted.out.bam",
        "bai": base_dir / "Aligned_sorted.out.bam.bai",
    }


def cdna_featurecounts_paths(sample):
    base_dir = cdna_star_dir(sample)
    return {
        "summary": base_dir / "featurecounts_summary.txt",
        "annotated": base_dir / "Aligned_sorted.out.bam.featureCounts.bam",
        "annotated_bai": base_dir / "Aligned_sorted.out.bam.featureCounts.bam.bai",
    }


def cdna_tagged_bam(sample):
    base_dir = cdna_star_dir(sample)
    return {
        "bam": base_dir / "Aligned_sorted_ftct_wTags.bam",
        "bai": base_dir / "Aligned_sorted_ftct_wTags.bam.bai",
    }


def cdna_counts(sample):
    return OUTPUT_DIR / "cdna" / sample / "counts.tsv.gz"


rule prepare_cdna_fastqs:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["cdna"]["r1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["cdna"]["r2"],
    output:
        r1=lambda wildcards: str(cdna_prepared_fastq(wildcards.sample, 1)),
        r2=lambda wildcards: str(cdna_prepared_fastq(wildcards.sample, 2)),
    params:
        trim=lambda wildcards: SAMPLES[wildcards.sample]["cdna"].get("trim", {}),
        threads=lambda wildcards: SAMPLES[wildcards.sample]["cdna"].get("trim", {}).get("threads", 4),
    run:
        from snakemake.shell import shell

        trim = params.trim
        out_dir = Path(cdna_prepared_fastq(wildcards.sample, 1)).parent
        out_dir.mkdir(parents=True, exist_ok=True)

        if trim.get("enabled", False):
            adapters = trim.get("adapters", {})
            adapters_r1 = " ".join(f"-a '{adapter}'" for adapter in adapters.get("r1", []))
            adapters_r2 = " ".join(f"-A '{adapter}'" for adapter in adapters.get("r2", []))
            minimum_length = trim.get("minimum_length", 20)
            quality_cutoff = trim.get("quality_cutoff", "20")
            extra_args = trim.get("extra_args", "")

            shell(
                "cutadapt -j {threads} {adapters_r1} {adapters_r2} "
                "-q {quality_cutoff} --minimum-length {minimum_length} {extra_args} "
                "-o {output.r1} -p {output.r2} {input.r1} {input.r2}"
            )
        else:
            for src, dest in [(input.r1, output.r1), (input.r2, output.r2)]:
                dest_path = Path(dest)
                dest_path.parent.mkdir(parents=True, exist_ok=True)
                if dest_path.exists() or dest_path.is_symlink():
                    dest_path.unlink()
                dest_path.symlink_to(Path(src).resolve())


rule parse_cdna_fastqs:
    input:
        r1=lambda wildcards: cdna_prepared_fastq(wildcards.sample, 1),
        r2=lambda wildcards: cdna_prepared_fastq(wildcards.sample, 2),
    output:
        fastq=lambda wildcards: str(cdna_barcoded_fastq(wildcards.sample)),
    params:
        genome_dir=lambda wildcards: config.get("split_pipe", {}).get("genome_dir"),
        chemistry=lambda wildcards: config.get("split_pipe", {}).get("chemistry", "v1"),
        extra=lambda wildcards: config.get("split_pipe", {}).get("extra_args", ""),
    shell:
        r"""
        set -euo pipefail
        tmpdir=$(mktemp -d)
        split-pipe --mode pre --chemistry {params.chemistry} \
            --genome_dir {params.genome_dir} \
            --fq1 {input.r1} --fq2 {input.r2} \
            --output_dir "$tmpdir" {params.extra}
        mkdir -p $(dirname {output.fastq})
        mv "$tmpdir"/process/barcode_head.fastq.gz {output.fastq}
        rm -rf "$tmpdir"
        """


rule add_cdna_indices:
    input:
        fastq=lambda wildcards: cdna_barcoded_fastq(wildcards.sample),
    output:
        fastq=lambda wildcards: str(cdna_indexed_fastq(wildcards.sample)),
    params:
        index=lambda wildcards: SAMPLES[wildcards.sample]["index"],
    shell:
        r"""
        mkdir -p $(dirname {output.fastq})
        python Scripts/combine_cell_bcs.py {input.fastq} {output.fastq} {params.index}
        """


rule star_align_cdna:
    input:
        fastq=lambda wildcards: cdna_indexed_fastq(wildcards.sample),
    output:
        bam=lambda wildcards: str(cdna_bam_paths(wildcards.sample)["bam"]),
        sorted=lambda wildcards: str(cdna_bam_paths(wildcards.sample)["sorted"]),
        bai=lambda wildcards: str(cdna_bam_paths(wildcards.sample)["bai"]),
    params:
        genome=lambda wildcards: config.get("star", {}).get("index"),
        extra=lambda wildcards: config.get("star", {}).get("extra_args", ""),
        prefix=lambda wildcards: str(cdna_star_dir(wildcards.sample)) + "/",
    threads=lambda wildcards: config.get("star", {}).get("threads", 8),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.prefix}
        STAR --runThreadN {threads} --genomeDir {params.genome} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} {params.extra}
        mv {params.prefix}Aligned.out.bam {output.bam}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        samtools index {output.sorted}
        """


rule feature_counts_cdna:
    input:
        bam=lambda wildcards: cdna_bam_paths(wildcards.sample)["sorted"],
    output:
        summary=lambda wildcards: str(cdna_featurecounts_paths(wildcards.sample)["summary"]),
        annotated=lambda wildcards: str(cdna_featurecounts_paths(wildcards.sample)["annotated"]),
        bai=lambda wildcards: str(cdna_featurecounts_paths(wildcards.sample)["annotated_bai"]),
    params:
        gtf=lambda wildcards: config.get("feature_counts", {}).get("gtf"),
        extra=lambda wildcards: config.get("feature_counts", {}).get("extra_args", "-t exon -g gene_id"),
    threads=lambda wildcards: config.get("feature_counts", {}).get("threads", 4),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.summary})
        featureCounts {params.extra} -a {params.gtf} -o {output.summary} -R BAM {input.bam}
        annotated_tmp="$(basename {input.bam}).featureCounts.bam"
        mv "$annotated_tmp" {output.annotated}
        samtools sort -@ {threads} -o {output.annotated}.tmp {output.annotated}
        mv {output.annotated}.tmp {output.annotated}
        samtools index {output.annotated}
        """


rule add_cdna_tags:
    input:
        bam=lambda wildcards: cdna_featurecounts_paths(wildcards.sample)["annotated"],
    output:
        bam=lambda wildcards: str(cdna_tagged_bam(wildcards.sample)["bam"]),
        bai=lambda wildcards: str(cdna_tagged_bam(wildcards.sample)["bai"]),
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        python Scripts/add_tags.py {input.bam} {output.bam}
        samtools index {output.bam}
        """


rule umi_tools_count_cdna:
    input:
        bam=lambda wildcards: cdna_tagged_bam(wildcards.sample)["bam"],
    output:
        counts=lambda wildcards: str(cdna_counts(wildcards.sample)),
    params:
        extra=lambda wildcards: config.get("umi_tools", {}).get("extra_args", "--per-gene --gene-tag=XT --assigned-status-tag=XS --extract-umi-method=tag --umi-tag=UB"),
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        umi_tools count {params.extra} -I {input.bam} -S {output.counts}
        """


rule combine_cdna_counts:
    input:
        lambda wildcards: [cdna_counts_path(sample) for sample in CDNA_SAMPLES],
    output:
        combined=lambda wildcards: combined_cdna_counts_path(),
    run:
        import gzip
        from pathlib import Path

        out_path = Path(output.combined)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        first = True
        with gzip.open(out_path, "wt") as out_handle:
            for path in input:
                with gzip.open(path, "rt") as in_handle:
                    if first:
                        for line in in_handle:
                            out_handle.write(line)
                        first = False
                    else:
                        next(in_handle)
                        for line in in_handle:
                            out_handle.write(line)
