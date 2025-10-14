ml load biology
ml load STAR

awk '
    /^>/ {
        if (keep) print seq;
        header=$0;
        match(header, /gene_symbol:([^ ]+)/, arr);
        gene=arr[1];
        match(header, /transcript_biotype:protein_coding/, type);
        biotype=(type[0] != "");
        keep = (gene ~ /^(AKIRIN1|ANXA1|APH1A|ARIH2|ATP2A3|BEX2|CTSL|DDIT4|EIF2A|HHEX|IARS2|KHDRBS1|MACF1|MFGE8|MRPL19|MTX1|MYC|PIM1|PRKAR2B|PRKCB|PRPF38A|RHAG|RPA2|SLC25A37|SS18L2|SSR3|XPO1)$/ && biotype);
        seq = "";
        if (keep) print header;
        next;
    }
    keep { seq = seq $0 }
    END { if (keep) print seq }
' ../../../ref_genomes/hg38_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa > STAR/selected_transcripts.fa

python Scripts/generate_gtf.py STAR/full_transcripts.fa STAR/full_transcripts.gtf

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR/full_transcripts --genomeFastaFiles STAR/full_transcripts.fa --sjdbGTFfile STAR/full_transcripts.gtf --sjdbOverhang 99
