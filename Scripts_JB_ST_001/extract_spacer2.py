import gzip
import sys

# Define flanking sequences for SPACER2
SPACER3_START = "TAGCTCTGAAAC"

TRNA_MAP = {
    "TGGAGGTACT": "Ala",
    "GGGCTCGTCC": "Pro",
    "AGGTTCCACC": "Gln"
}

def extract_between(seq, start, trna_map):
    """Extracts the sequence between two defined flanking sequences."""
    start_idx = seq.find(start)
    if start_idx == -1:
        return "NOSEQ"

    min_end_idx = None
    for trna_seq in trna_map.keys():
        end_idx = seq.find(trna_seq, start_idx + len(start))  # Ensure end comes after start
        if end_idx != -1:
            if min_end_idx is None or end_idx < min_end_idx:
                min_end_idx = end_idx

    if min_end_idx is not None:
        return seq[start_idx + len(start) : min_end_idx]  # Extract in-between sequence
    return "NOSEQ"

def find_trna_id(seq):
    """Search for a known tRNA sequence in the read and return the corresponding tRNA_id."""
    for trna_seq, trna_id in TRNA_MAP.items():
        if trna_seq in seq:
            return trna_id
    return "UNKNOWN"  # If no match is found

def modify_fastq_headers(barcoded_r1_file, output_r1):
    """Modifies R1 headers to include extracted SPACER3 sequence."""
    count_extracted = 0  # Counter for successful extractions

    with gzip.open(barcoded_r1_file, "rt") as r1, gzip.open(output_r1, "wt") as out_r1:
        while True:
            # Read 4 lines per read (FASTQ format)
            r1_lines = [r1.readline().strip() for _ in range(4)]

            if not r1_lines[0]:
                break  # Stop at EOF

            # Extract the new spacer (SPACER3) from the R1 sequence
            spacer3_seq = extract_between(r1_lines[1], SPACER3_START, TRNA_MAP)
            trna_id = find_trna_id(r1_lines[1])

            # Check if extraction was successful
            if spacer3_seq != "NOSEQ":
                count_extracted += 1  # Increment count

            # Modify R1 header to include SPACER3
            new_header = r1_lines[0] + f":SPACER2={spacer3_seq}:tRNA={trna_id}"
            r1_lines[0] = new_header  # Update header

            # Write modified R1 to output
            out_r1.write("\n".join(r1_lines) + "\n")

    # Print total number of extracted SPACER3 sequences
    print(f"Total successfully extracted SPACER2 sequences: {count_extracted}")

# Command-line execution
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_spacer2.py input_barcoded_R1.fastq.gz output_R1.fastq.gz")
        sys.exit(1)

    modify_fastq_headers(sys.argv[1], sys.argv[2])

