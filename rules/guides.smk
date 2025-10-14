from pathlib import Path


def guide_mode(sample):
    return SAMPLES[sample]["guides"].get("mode", "dual").lower()


def _guide_prepare_dir(sample):
    return OUTPUT_DIR / "intermediate" / sample / "guides"


def guide_prepared_fastq(sample, read):
    return _guide_prepare_dir(sample) / f"{sample}_guides_R{read}_prepared.fastq.gz"


def guide_barcoded_fastq(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_barcode_head.fastq.gz"


def guide_indexed_fastq(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_barcode_head_wIndex.fastq.gz"


def guide_spacer1_fastq(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_barcode_head_spacer1.fastq.gz"


def guide_spacers_fastq(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_barcode_head_wSpacers.fastq.gz"


def guide_assignments(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_assignments.tsv"


def guide_summary(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_summary.txt"


def guide_read_counts(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_read_counts.tsv"


def guide_collapsed(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_assignments_collapsed.tsv"


def guide_filtered(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_assignments_filtered.tsv"


def guide_counts(sample):
    return OUTPUT_DIR / "guides" / sample / "guide_counts.tsv"


def guide_umi_log(sample):
    return OUTPUT_DIR / "guides" / sample / "umi_corrections.tsv"


def guide_spacer_log(sample):
    return OUTPUT_DIR / "guides" / sample / "spacer_corrections.tsv"


def _require_dual(sample):
    if sample not in DUAL_GUIDE_SAMPLES:
        raise ValueError(f"Sample {sample} is not configured as dual-guide")


def _require_single(sample):
    if sample not in SINGLE_GUIDE_SAMPLES:
        raise ValueError(f"Sample {sample} is not configured as single-guide")


def dual_indexed_fastq(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_indexed_fastq(wildcards.sample))


def dual_spacer1_fastq(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_spacer1_fastq(wildcards.sample))


def dual_spacers_fastq(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_spacers_fastq(wildcards.sample))


def dual_assignments_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_assignments(wildcards.sample))


def dual_summary_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_summary(wildcards.sample))


def dual_read_counts_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_read_counts(wildcards.sample))


def dual_collapsed_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_collapsed(wildcards.sample))


def dual_filtered_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_filtered(wildcards.sample))


def dual_counts_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_counts(wildcards.sample))


def dual_umi_log_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_umi_log(wildcards.sample))


def dual_spacer_log_path(wildcards):
    _require_dual(wildcards.sample)
    return str(guide_spacer_log(wildcards.sample))


def single_indexed_fastq(wildcards):
    _require_single(wildcards.sample)
    return str(guide_indexed_fastq(wildcards.sample))


def single_assignments_path(wildcards):
    _require_single(wildcards.sample)
    return str(guide_assignments(wildcards.sample))


def single_collapsed_path(wildcards):
    _require_single(wildcards.sample)
    return str(guide_collapsed(wildcards.sample))


def single_filtered_path(wildcards):
    _require_single(wildcards.sample)
    return str(guide_filtered(wildcards.sample))


def single_counts_path(wildcards):
    _require_single(wildcards.sample)
    return str(guide_counts(wildcards.sample))


if GUIDE_SAMPLES:

    rule prepare_guide_fastqs:
        input:
            r1=lambda wildcards: SAMPLES[wildcards.sample]["guides"]["r1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["guides"]["r2"],
        output:
            r1=lambda wildcards: str(guide_prepared_fastq(wildcards.sample, 1)),
            r2=lambda wildcards: str(guide_prepared_fastq(wildcards.sample, 2)),
        params:
            config=lambda wildcards: SAMPLES[wildcards.sample]["guides"],
        run:
            from snakemake.shell import shell
            guides_cfg = params.config
            out_dir = Path(guide_prepared_fastq(wildcards.sample, 1)).parent
            out_dir.mkdir(parents=True, exist_ok=True)

            if guides_cfg.get("spacer_to_r1", False):
                shell(
                    "python Scripts/spacer_to_R1header.py {input.r1} {input.r2} {output.r1} {output.r2}"
                )
            else:
                for src, dest in [(input.r1, output.r1), (input.r2, output.r2)]:
                    dest_path = Path(dest)
                    dest_path.parent.mkdir(parents=True, exist_ok=True)
                    if dest_path.exists() or dest_path.is_symlink():
                        dest_path.unlink()
                    dest_path.symlink_to(Path(src).resolve())


    rule parse_guide_fastqs:
        input:
            r1=lambda wildcards: guide_prepared_fastq(wildcards.sample, 1),
            r2=lambda wildcards: guide_prepared_fastq(wildcards.sample, 2),
        output:
            fastq=lambda wildcards: str(guide_barcoded_fastq(wildcards.sample)),
        params:
            genome_dir=lambda wildcards: config.get("split_pipe", {}).get("genome_dir"),
            chemistry=lambda wildcards: SAMPLES[wildcards.sample]["guides"].get("chemistry", config.get("split_pipe", {}).get("chemistry", "v1")),
            extra=lambda wildcards: SAMPLES[wildcards.sample]["guides"].get("split_pipe_extra", config.get("split_pipe", {}).get("extra_args", "")),
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


    rule add_guide_indices:
        input:
            fastq=lambda wildcards: guide_barcoded_fastq(wildcards.sample),
        output:
            fastq=lambda wildcards: str(guide_indexed_fastq(wildcards.sample)),
        params:
            index=lambda wildcards: SAMPLES[wildcards.sample]["index"],
        shell:
            r"""
            mkdir -p $(dirname {output.fastq})
            python Scripts/combine_cell_bcs.py {input.fastq} {output.fastq} {params.index}
            """


if DUAL_GUIDE_SAMPLES:

    rule annotate_spacer1:
        input:
            fastq=dual_indexed_fastq,
        output:
            fastq=dual_spacer1_fastq,
        shell:
            r"""
            mkdir -p $(dirname {output.fastq})
            bash Scripts/spacer_to_header.sh {input.fastq} {output.fastq}
            """


    rule extract_spacer2:
        input:
            fastq=dual_spacer1_fastq,
        output:
            fastq=dual_spacers_fastq,
        shell:
            r"""
            python Scripts/extract_spacer2.py {input.fastq} {output.fastq}
            """


    rule extract_guides_dual:
        input:
            fastq=dual_spacers_fastq,
        output:
            assignments=dual_assignments_path,
            summary=dual_summary_path,
            counts=dual_read_counts_path,
        shell:
            r"""
            mkdir -p $(dirname {output.assignments})
            python Scripts/extract_guides.py {input.fastq} {output.assignments} {output.summary} {output.counts}
            """


    rule collapse_dual_guides:
        input:
            counts=dual_read_counts_path,
        output:
            collapsed=dual_collapsed_path,
            umi_log=dual_umi_log_path,
            spacer_log=dual_spacer_log_path,
        shell:
            r"""
            mkdir -p $(dirname {output.collapsed})
            python Scripts/collapse_umis.py {input.counts} {output.collapsed} {output.umi_log} {output.spacer_log}
            """


    rule filter_dual_guides:
        input:
            collapsed=dual_collapsed_path,
        output:
            filtered=dual_filtered_path,
        params:
            threshold=lambda wildcards: str(SAMPLES[wildcards.sample]["guides"].get("tpt_threshold", config.get("guide_processing", {}).get("tpt_threshold", 0.25))),
            histogram=lambda wildcards: histogram_path(wildcards.sample),
        shell:
            r"""
            mkdir -p $(dirname {output.filtered})
            python Scripts/theseus_impl.py {input.collapsed} {output.filtered} {params.threshold} {params.histogram}
            """


    rule summarize_dual_guides:
        input:
            filtered=dual_filtered_path,
        output:
            counts=dual_counts_path,
        shell:
            r"""
            python Scripts/summarize_dual_guides.py {input.filtered} {output.counts}
            """


if SINGLE_GUIDE_SAMPLES:

    rule process_single_guides:
        input:
            fastq=single_indexed_fastq,
        output:
            assignments=single_assignments_path,
            aggregated=single_collapsed_path,
        params:
            reference=lambda wildcards: SAMPLES[wildcards.sample]["guides"].get("reference", config.get("guide_processing", {}).get("single_reference")),
        shell:
            r"""
            mkdir -p $(dirname {output.assignments})
            python Scripts/process_single_guides.py {input.fastq} {params.reference} {output.assignments} {output.aggregated}
            """

    rule filter_single_guides:
        input:
            aggregated=single_collapsed_path,
        output:
            filtered=single_filtered_path,
        params:
            threshold=lambda wildcards: str(SAMPLES[wildcards.sample]["guides"].get("tpt_threshold", config.get("guide_processing", {}).get("tpt_threshold", 0.25))),
            histogram=lambda wildcards: histogram_path(wildcards.sample),
        shell:
            r"""
            mkdir -p $(dirname {output.filtered})
            python Scripts/theseus_impl.py {input.aggregated} {output.filtered} {params.threshold} {params.histogram}
            """

    rule summarize_single_guides:
        input:
            filtered=single_filtered_path,
        output:
            counts=single_counts_path,
        shell:
            r"""
            python Scripts/summarize_single_guides.py {input.filtered} {output.counts}
            """


if GUIDE_SAMPLES:

    rule combine_guide_counts:
        input:
            lambda wildcards: [guide_counts_path(sample) for sample in GUIDE_SAMPLES],
        output:
            combined=lambda wildcards: combined_guide_counts_path(),
        run:
            import csv
            from pathlib import Path

            out_path = Path(output.combined)
            out_path.parent.mkdir(parents=True, exist_ok=True)

            fieldnames = []
            rows = []
            for path in input:
                with open(path, "r", newline="") as in_handle:
                    reader = csv.DictReader(in_handle, delimiter="\t")
                    if reader.fieldnames is None:
                        continue
                    for name in reader.fieldnames:
                        if name not in fieldnames:
                            fieldnames.append(name)
                    rows.extend(list(reader))

            with out_path.open("w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow({name: row.get(name, "") for name in fieldnames})
