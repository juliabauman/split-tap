# Rules for combining final outputs

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
