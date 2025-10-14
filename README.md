# Split-TAP Snakemake Workflow

This repository contains a Snakemake workflow that standardises the Split-TAP processing steps for both cDNA and guide libraries. The pipeline loads its configuration from `config/config.yaml` and exposes common helper functions from the top-level `Snakefile` so that rule modules under `rules/` can build per-sample targets automatically.

## Repository layout

- `Snakefile` – entry point that wires the configuration to the cDNA and guide rule modules and defines the default targets.  
- `rules/cdna.smk` – rules for preparing, parsing, aligning, and counting cDNA reads.  
- `rules/guides.smk` – rules for processing single- and dual-guide libraries, merging count tables, and generating histograms.  
- `config/config.yaml` – tool defaults (Split-Pipe, STAR, featureCounts, UMI-tools) and per-sample metadata, including which analyses should run.  
- `Scripts/` – helper scripts invoked by the Snakemake rules (for example `combine_cell_bcs.py`, `process_single_guides.py`, and `summarize_dual_guides.py`).

Refer to each file for additional details about paths, parameters, and expected inputs.

## Running the workflow

1. Update `config/config.yaml` with the FASTQ locations, barcode index, genome references, and other parameters for your experiment. Samples can independently enable or disable cDNA or guide branches, and guide samples can specify `mode: dual` or `mode: single` to pick the correct rules.
2. Create the output directory specified by `paths.output_dir` if it does not already exist.
3. Launch Snakemake from the repository root, for example:
   ```bash
   snakemake --cores 8
   ```
   Adjust the `--cores` value to match your environment. Snakemake will materialise the per-sample and combined outputs declared in the `rule all` inputs.

The workflow expects the command-line tools referenced in the rules (Cutadapt, split-pipe, STAR, samtools, featureCounts, umi_tools, and Python) to be available on your `PATH`.

## Pushing changes to GitHub

This environment cannot authenticate to your GitHub remote, so pushing must be done from your machine. After cloning the repository locally and pulling the latest commits, run:

```bash
git push origin <branch>
```

Replace `<branch>` with `main` or whichever branch you are using. Set up the remote with `git remote add origin <url>` if it has not been configured yet.
