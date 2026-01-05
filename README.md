# Perturb-seq Analysis Pipeline

A Snakemake pipeline for processing Perturb-seq data with support for both single and paired guide modes.

## Overview

This pipeline processes Perturb-seq data from raw FASTQ files through to final count matrices, supporting:
- **Single guide mode**: One guide per construct (e.g., promoter targeting)
- **Paired guide mode**: Two guides per construct (e.g., enhancer pairs)

The pipeline handles both cDNA (for gene expression) and CROP (for guide assignment) data in parallel.

## Directory Structure

```
.
├── Snakefile              # Main pipeline file
├── config.yaml            # Configuration file
├── environment.yml        # Conda environment specification
├── rules/                 # Snakemake rule files
│   ├── cdna_processing.smk
│   ├── crop_processing.smk
│   ├── star_alignment.smk
│   └── final_outputs.smk
├── Scripts/               # Analysis scripts (symlinked from previous analyses)
├── i7_barcodes.csv        # Sample barcoding information
├── all_spacer_pairs_rc.csv      # Guide pair library (paired mode)
└── targetingGuides_Final.csv    # Guide library (single mode)
```

## Setup

### 1. Environment Setup

This pipeline requires a conda environment with all necessary dependencies. Create the environment from the included specification:

```bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate split-tap
```

The environment includes:
- **split-pipe** (v1.1.1) - SPLiT-seq barcode parsing
- **cutadapt** (v4.6) - Adapter trimming
- **STAR** - RNA-seq alignment (via subread package)
- **samtools** - BAM file manipulation
- **featureCounts** (subread v2.0.8) - Gene counting
- **umi_tools** (v1.1.6) - UMI deduplication
- **pigz** (v2.8) - Parallel gzip
- **Python 3.10** with pandas, pysam, biopython
- **scanpy** and analysis packages

**Note:** You must activate this environment before running the pipeline.

### 2. Input Data Structure

Organize your FASTQ files in the following structure:
```
data_dir/
  ├── WELL_1/    # cDNA files (e.g., A1_1/)
  │   ├── *_R1.fastq.gz
  │   └── *_R2.fastq.gz
  └── WELL_2/    # CROP files (e.g., A1_2/)
      ├── *_R1.fastq.gz
      └── *_R2.fastq.gz
```

Where `WELL` matches entries in your `config.yaml` (e.g., A1, A2, etc.)

### 3. Configuration Files

#### config.yaml
Edit the main configuration file to specify:
- `guide_mode`: "single" or "paired"
- `data_dir`: Path to raw FASTQ files
- `wells`: List of wells to process (e.g., ["A1", "A3", "A4"])
- `genome_dir`: Path to reference genome
- `star_genome_dir`: Path to STAR index
- `star_gtf`: Path to GTF file
- Thread and memory settings

#### i7_barcodes.csv
Format:
```csv
well,sample_name,i7_index
A1,SPRITE_PCR_P7_1,AAGTAGAG
A3,SPRITE_PCR_P7_3,TGTTCCGA
```

#### Guide Library Files

**For paired mode** - `all_spacer_pairs_rc.csv`:
```csv
Guide1-Guide2_Name,Seq1,Seq2,Seq2_rc,Seq1-Seq2_rc
enh1_guide1-enh2_guide1,GGGCAACTGCTTATCTAACAT,GCCAAACCACAGGCTCCCTCA,TGAGGGAGCCTGTGGTTTGGC,GGGCAACTGCTTATCTAACATTGAGGGAGCCTGTGGTTTGGC
```

**For single mode** - `targetingGuides_Final.csv`:
```csv
type,chromosome,start,end,guide_name,score,strand,guide_seq,guide_seq_wG,PAM,position
promoter_targeting,chr7,107044723,107044742,GENE1_guide1,0,+,ACACGGAGCAGACGCGCGCC,GACACGGAGCAGACGCGCGCC,NNN,chr7:107044723-107044742
```

## Running the Pipeline

**Important:** Always activate the conda environment first:
```bash
conda activate split-tap
```

### Dry Run (Test)
```bash
snakemake -n
```

### Local Execution
```bash
snakemake --cores 8
```

### SLURM Cluster Execution
```bash
snakemake --cluster "sbatch --partition=larsms --cpus-per-task={threads} --mem={resources.mem_mb}M --time=12:00:00" \
          --jobs 20 \
          --cores 100
```

### With Snakemake Profile (Recommended)
Create a profile for your cluster and run:
```bash
snakemake --profile slurm --jobs 20
```

## Pipeline Workflow

### Common Steps (Both Modes)
1. **cDNA Processing**:
   - Parse barcoding (split-pipe)
   - Cutadapt adapter/polyA trimming
   - Add i7 index to cell barcodes
   - STAR alignment
   - Add CB/UB tags
   - featureCounts
   - UMI counting with umi_tools

### Paired Guide Mode
2. **CROP Processing**:
   - Extract spacer from R2 to R1 header
   - Parse barcoding
   - Add i7 index
   - Extract first spacer to header
   - Extract second spacer to header
   - Combine all wells
   - Extract guide information to table
   - Split into chunks for parallel processing
   - Collapse UMIs (parallel)
   - Chimera correction (theseus)
   - Filter by chimera results
   - Assign guide pair IDs
   - Assign tRNA for covariate

### Single Guide Mode
2. **CROP Processing**:
   - Parse barcoding
   - Add i7 index
   - Combine all wells
   - Extract single guides
   - Chimera correction (theseus)
   - Count guides

## Output Files

### Final Outputs (Data/)
- `combined_cDNA_counts.tsv` - Combined gene expression counts across all wells
- **Paired mode**:
  - `CROP_counts.txt` - Guide pair counts per cell
  - `tRNA_assignments.txt` - tRNA assignments for covariates
- **Single mode**:
  - `guide_counts.txt` - Single guide counts per cell

### Intermediate Files
- `barcoded_fastqs/` - Barcoded FASTQ files per well
- `star_outputs/` - STAR alignment results per well
- `featureCounts/` - featureCounts outputs per well
- `final_outputs/` - Per-well cDNA counts

## Key Features

1. **Flexible Input**: Supports variable numbers of wells and samples
2. **Efficient Processing**: Parallelizes across wells and uses checkpoints for dynamic chunking
3. **Resource Management**: Configurable threads and memory per rule
4. **Mode Switching**: Single configuration file controls single vs paired guide workflows
5. **Adapter Trimming**: Always applies cutadapt to cDNA reads for quality
6. **Chimera Correction**: Applies theseus algorithm to filter chimeric molecules

## Troubleshooting

### Missing Input Files
Ensure your data directory structure matches the expected format (WELL_1 and WELL_2 folders)

### Module Conflicts
The pipeline uses various modules. If you encounter conflicts:
- Purge modules before running: `module purge`
- Let Snakemake load required modules per rule

### Memory Issues
Increase memory in `config.yaml`:
- `resources.parse_mem`: For split-pipe
- `resources.star_mem`: For STAR alignment
- `resources.featurecounts_mem`: For featureCounts

### Thread Allocation
Adjust thread counts in `config.yaml` based on available resources

## Notes

- The pipeline uses temporary files extensively to save disk space
- CROP UMI collapse is parallelized into 58 chunks by default
- cDNA cutadapt trimming is always enabled (removes adapters and polyA tails)
- All required software is included in the `environment.yml` conda environment

## Contact

For questions or issues, consult the original analysis scripts in:
- `JB_ST_001/screen_reprocess/Scripts/` (single guide mode)
- `JB_ST_002/screen/Scripts/` (paired guide mode)
