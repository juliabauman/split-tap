# Rules for combining final outputs and SCEPTRE preprocessing

rule combine_cdna_counts:
    """Combine cDNA counts from all wells"""
    input:
        expand(f"{OUTPUT_DIR}/final_outputs/{{well}}_cDNA_counts.tsv", well=CDNA_WELLS)
    output:
        f"{OUTPUT_DIR}/combined_cDNA_counts.tsv"
    run:
        if len(input) > 0:
            shell("head -n 1 {input[0]} > {output}")
            shell("tail -n +2 -q {input} >> {output}")
        else:
            shell("touch {output}")


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
            # Columns: cell, Guide1-Guide2_Name, UMI_Count
            guides_df.rename(columns={'Guide1-Guide2_Name':'gene', 'UMI_Count':'count'}, inplace=True)
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
        "../Scripts/create_knee_plot.R"
