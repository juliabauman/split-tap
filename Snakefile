configfile: "config/config.yaml"

from pathlib import Path

SAMPLES = config.get("samples", {})
OUTPUT_DIR = Path(config.get("paths", {}).get("output_dir", "Data"))

CDNA_SAMPLES = [sample for sample, cfg in SAMPLES.items() if cfg.get("cdna") and cfg.get("cdna", {}).get("enabled", True)]
GUIDE_SAMPLES = [sample for sample, cfg in SAMPLES.items() if cfg.get("guides") and cfg.get("guides", {}).get("enabled", True)]
DUAL_GUIDE_SAMPLES = [sample for sample in GUIDE_SAMPLES if SAMPLES[sample]["guides"].get("mode", "dual").lower() == "dual"]
SINGLE_GUIDE_SAMPLES = [sample for sample in GUIDE_SAMPLES if SAMPLES[sample]["guides"].get("mode", "dual").lower() != "dual"]


def cdna_counts_path(sample):
    return str(OUTPUT_DIR / "cdna" / sample / "counts.tsv.gz")


def guide_counts_path(sample):
    return str(OUTPUT_DIR / "guides" / sample / "guide_counts.tsv")


def combined_cdna_counts_path():
    return str(OUTPUT_DIR / "cdna" / "combined_counts.tsv.gz")


def combined_guide_counts_path():
    return str(OUTPUT_DIR / "guides" / "combined_counts.tsv")


def histogram_path(sample):
    hist_dir = Path(config.get("guide_processing", {}).get("histogram_dir", OUTPUT_DIR / "guides" / "histograms"))
    return str(Path(hist_dir) / f"{sample}_tpt_histogram.tsv")


workflow.globals.update(
    {
        "OUTPUT_DIR": OUTPUT_DIR,
        "SAMPLES": SAMPLES,
        "CDNA_SAMPLES": CDNA_SAMPLES,
        "GUIDE_SAMPLES": GUIDE_SAMPLES,
        "DUAL_GUIDE_SAMPLES": DUAL_GUIDE_SAMPLES,
        "SINGLE_GUIDE_SAMPLES": SINGLE_GUIDE_SAMPLES,
        "cdna_counts_path": cdna_counts_path,
        "guide_counts_path": guide_counts_path,
        "combined_cdna_counts_path": combined_cdna_counts_path,
        "combined_guide_counts_path": combined_guide_counts_path,
        "histogram_path": histogram_path,
    }
)

include: "rules/cdna.smk"
include: "rules/guides.smk"


rule all:
    input:
        ([combined_cdna_counts_path()] if CDNA_SAMPLES else [])
        + [cdna_counts_path(sample) for sample in CDNA_SAMPLES]
        + ([combined_guide_counts_path()] if GUIDE_SAMPLES else [])
        + [guide_counts_path(sample) for sample in GUIDE_SAMPLES]
