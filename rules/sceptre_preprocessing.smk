# Rules for SCEPTRE preprocessing
# Based on make_sceptre_files.R preprocessing section

rule combine_guide_cdna_counts:
    """Combine guide and cDNA counts into single file for SCEPTRE"""
    input:
        cdna = f"{OUTPUT_DIR}/combined_cDNA_counts.tsv",
        guides = f"{OUTPUT_DIR}/CROP_counts.txt" if GUIDE_MODE == "paired" else f"{OUTPUT_DIR}/guide_counts.txt"
    output:
        f"{OUTPUT_DIR}/guide_cDNA_counts.tsv"
    run:
        import pandas as pd

        # Read cDNA data (columns: gene, cell, count)
        cdna_df = pd.read_csv(input.cdna, sep='\t', header=0)
        # Rename to match guide format
        cdna_df = cdna_df[['cell', 'gene', 'count']]

        # Read guide data
        guides_df = pd.read_csv(input.guides, sep='\t')

        # Standardize guide column names based on mode
        if GUIDE_MODE == "paired":
            # Columns: cell, Guide1_Guide2_Name, UMI_Count
            guides_df.rename(columns={'Guide1_Guide2_Name':'gene', 'UMI_Count':'count'}, inplace=True)
        else:
            # Single mode - adjust column names as needed
            # Assuming columns: Cell_Barcode, Guide_Name, UMI_Count
            if 'gRNA_name' in guides_df.columns:
                guides_df.rename(columns={'gRNA_name': 'gene'}, inplace=True)
            elif 'guide_name' in guides_df.columns:
                guides_df.rename(columns={'guide_name': 'gene'}, inplace=True)

        guides_df = guides_df[['cell', 'gene', 'count']]

        # Combine both datasets
        combined = pd.concat([guides_df, cdna_df], ignore_index=True)

        # Write to output
        combined.to_csv(output[0], sep='\t', index=False)


rule create_knee_plot:
    """Generate knee plot to determine cell filtering threshold"""
    input:
        f"{OUTPUT_DIR}/combined_cDNA_counts.tsv"
    output:
        plot = f"{OUTPUT_DIR}/plots/cell_knee_plot.pdf",
        summary = f"{OUTPUT_DIR}/cell_umi_summary.tsv"
    params:
        test_thresh = config.get("sceptre", {}).get("knee_plot_threshold", 700000)
    script:
        "../Scripts//create_knee_plot.R"


rule filter_cells_by_knee:
    """Filter cells based on knee threshold and create filtered datasets"""
    input:
        guide_cdna = f"{OUTPUT_DIR}/guide_cDNA_counts.tsv",
        guides = f"{OUTPUT_DIR}/CROP_counts.txt" if GUIDE_MODE == "paired" else f"{OUTPUT_DIR}/guide_counts.txt",
        cdna = f"{OUTPUT_DIR}/combined_cDNA_counts.tsv",
        umi_summary = f"{OUTPUT_DIR}/cell_umi_summary.tsv"
    output:
        filtered_guide_cdna = f"{OUTPUT_DIR}/guide_cDNA_counts_filtered.tsv",
        filtered_cells = f"{OUTPUT_DIR}/filtered_cell_list.txt",
        umi_histogram = f"{OUTPUT_DIR}/plots/kept_cells_umi_histogram.pdf",
        stats = f"{OUTPUT_DIR}/cell_filtering_stats.txt"
    params:
        cell_thresh = config.get("sceptre", {}).get("cell_threshold", 700000),
        max_umi_count = config.get("sceptre", {}).get("max_umi_per_cell", 300)
    script:
        "../Scripts/filter_cells_by_knee.R"


rule create_sceptre_inputs:
    """Create SCEPTRE input files (matrix.mtx, features.tsv, barcodes.tsv)"""
    input:
        guide_cdna = f"{OUTPUT_DIR}/guide_cDNA_counts_filtered.tsv",
        filtered_cells = f"{OUTPUT_DIR}/filtered_cell_list.txt"
    output:
        matrix = f"{OUTPUT_DIR}/sceptre_files/matrix.mtx",
        features = f"{OUTPUT_DIR}/sceptre_files/features.tsv",
        barcodes = f"{OUTPUT_DIR}/sceptre_files/barcodes.tsv"
    params:
        gene_list = ['AKIRIN1', 'ANXA1', 'APH1A', 'ARIH2', 'ARL4A', 'ATP2A3', 'BEX2', 'CTSL', 'DDIT4',
                     'EIF2A', 'HHEX', 'IARS2', 'KHDRBS1', 'MYC', 'PIM1', 'MACF1', 'MFGE8', 'MRPL19',
                     'MTX1', 'NET1', 'PRKAR2B', 'PRKCB', 'PRPF38A', 'RHAG', 'RPA2', 'SLC25A37',
                     'SS18L2', 'SSR3', 'XPO1']
    envmodules: "R/4.4.2"
    script:
        "../Scripts/create_sceptre_inputs.R"
