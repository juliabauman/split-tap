configfile: "config/config.yaml"

rule all:
    input:
        "../Data/cdna_counts.tsv.gz"

include: "rules/preprocess.smk"
include: "rules/count.smk"
