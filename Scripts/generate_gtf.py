import sys

def parse_fasta(fasta_file):
    """Reads a FASTA file and returns a dictionary of transcript lengths."""
    transcript_lengths = {}
    with open(fasta_file, "r") as f:
        transcript_id = None
        sequence = []
        for line in f:
            if line.startswith(">"):
                if transcript_id and sequence:
                    transcript_lengths[transcript_id] = len("".join(sequence))
                transcript_id = line.split()[0][1:]  # Extract transcript ID (remove ">")
                sequence = []
            else:
                sequence.append(line.strip())

        if transcript_id and sequence:
            transcript_lengths[transcript_id] = len("".join(sequence))

    return transcript_lengths

def generate_gtf(transcript_lengths, output_gtf):
    """Creates a GTF file with transcript-relative coordinates."""
    with open(output_gtf, "w") as gtf_file:
        for transcript_id, length in transcript_lengths.items():
            # Assume all transcripts are on the + strand (modify if necessary)
            strand = "+"
            gene_id = f"{transcript_id}_GENE"
            gene_name = f"{transcript_id}_GENE"

            # GTF formatted entries
            gene_entry = f"{transcript_id}\tcustom\tgene\t1\t{length}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_name}\";\n"
            transcript_entry = f"{transcript_id}\tcustom\ttranscript\t1\t{length}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\"; gene_name \"{gene_name}\";\n"
            exon_entry = f"{transcript_id}\tcustom\texon\t1\t{length}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\"; gene_name \"{gene_name}\"; exon_number \"1\";\n"

            # Write to file
            gtf_file.writelines([gene_entry, transcript_entry, exon_entry])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_transcript_gtf.py <input_fasta> <output_gtf>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_gtf = sys.argv[2]

    transcript_lengths = parse_fasta(input_fasta)
    generate_gtf(transcript_lengths, output_gtf)

    print(f"GTF file generated: {output_gtf}")

