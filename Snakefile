import os
import glob
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Extract configuration parameters
GUIDE_MODE = config["guide_mode"]
DATA_DIR = config["data_dir"]
WELLS = config["wells"]
OUTPUT_DIR = config["output_dir"]
SCRIPTS_DIR = config["scripts_dir"]
GENOME_DIR = config["genome_dir"]
STAR_GENOME_DIR = config["star_genome_dir"]
STAR_GTF = config["star_gtf"]

# Validate guide mode
if GUIDE_MODE not in ["single", "paired"]:
    raise ValueError(f"guide_mode must be 'single' or 'paired', got '{GUIDE_MODE}'")

# Helper function to find input folders
def get_well_folders():
    """Find all well folders matching the pattern WELL_1 and WELL_2"""
    folders = {"cdna": {}, "crop": {}}
    for well in WELLS:
        # Find cDNA folders (ending in _1)
        cdna_pattern = os.path.join(DATA_DIR, f"{well}_1")
        if os.path.isdir(cdna_pattern):
            folders["cdna"][well] = cdna_pattern

        # Find CROP folders (ending in _2)
        crop_pattern = os.path.join(DATA_DIR, f"{well}_2")
        if os.path.isdir(crop_pattern):
            folders["crop"][well] = crop_pattern

    return folders

WELL_FOLDERS = get_well_folders()
CDNA_WELLS = list(WELL_FOLDERS["cdna"].keys())
CROP_WELLS = list(WELL_FOLDERS["crop"].keys())

# Define final outputs based on guide mode
def get_final_outputs():
    outputs = []

    # Combined cDNA counts (always generated)
    outputs.append(f"{OUTPUT_DIR}/combined_cDNA_counts.tsv")

    # Guide counts depend on mode
    if GUIDE_MODE == "paired":
        outputs.append(f"{OUTPUT_DIR}/CROP_counts.txt")
        outputs.append(f"{OUTPUT_DIR}/tRNA_assignments.txt")
    else:  # single mode
        outputs.append(f"{OUTPUT_DIR}/guide_counts.txt")

    # SCEPTRE preprocessing outputs
    outputs.append(f"{OUTPUT_DIR}/guide_cDNA_counts.tsv")
    outputs.append(f"{OUTPUT_DIR}/plots/cell_knee_plot.pdf")
    outputs.append(f"{OUTPUT_DIR}/guide_cDNA_counts_filtered.tsv")
    outputs.append(f"{OUTPUT_DIR}/plots/kept_cells_umi_histogram.pdf")
    outputs.append(f"{OUTPUT_DIR}/sceptre_files/matrix.mtx")
    outputs.append(f"{OUTPUT_DIR}/sceptre_files/features.tsv")
    outputs.append(f"{OUTPUT_DIR}/sceptre_files/barcodes.tsv")

    return outputs

# Main rule - runs full pipeline including SCEPTRE preprocessing
rule all:
    input:
        get_final_outputs()

# Inspection rule - generates only the diagnostic plots for manual threshold selection
rule inspect_plots:
    input:
        f"{OUTPUT_DIR}/guide_cDNA_counts.tsv",
        f"{OUTPUT_DIR}/plots/cell_knee_plot.pdf",
        f"{OUTPUT_DIR}/cell_umi_summary.tsv"

# Include rule files
include: "rules/cdna_processing.smk"
include: "rules/crop_processing.smk"
include: "rules/star_alignment.smk"
include: "rules/final_outputs.smk"
include: "rules/sceptre_preprocessing.smk"
