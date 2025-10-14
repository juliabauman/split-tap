import gzip
import sys
import re

# Define the flanking sequences
START_SEQ = "CCAATGCA"
END_SEQ = "GTTTGAGAGCTAAGCA"
LINKER_SEQ="LINKER"

def extract_between(seq, start, end):
    """Extracts the sequence between two defined flanking sequences."""
    start_idx = seq.find(start)
    end_idx = seq.find(end, start_idx + len(start))  # Ensure end comes after start

    if start_idx != -1 and end_idx != -1:
        return seq[start_idx + len(start) : end_idx]  # Extract in-between sequence
    return None

def modify_fastq_headers(r1_file, r2_file, output_r1, output_r2):
    """Modifies R1 headers to include extracted sequence from R2."""
    count_extracted = 0 
    with gzip.open(r1_file, "rt") as r1, gzip.open(r2_file, "rt") as r2, gzip.open(output_r1, "wt") as out_r1, gzip.open(output_r2, "wt") as out_r2:
        while True:
            # Read 4 lines per read (FASTQ format)
            r1_lines = [r1.readline().strip() for _ in range(4)]
            r2_lines = [r2.readline().strip() for _ in range(4)]
            
            # Stop at EOF
            if not r1_lines[0]:
                break
            
            # Extract the target sequence from R2
            spacer = extract_between(r2_lines[1], START_SEQ, END_SEQ)

            if spacer:
                count_extracted += 1
                r1_lines[1] += LINKER_SEQ + spacer
                r1_lines[3] += "I" * (len(LINKER_SEQ) + len(spacer))

            # Write modified R1 to output
            out_r1.write("\n".join(r1_lines) + "\n")
            out_r2.write("\n".join(r2_lines) + "\n")

    print(f"Total successfully extracted sequences: {count_extracted}")

# Command-line execution
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python spacer_to_R1.py input_R1.fastq.gz input_R2.fastq.gz output_R1.fastq.gz output_R2.fastq.gz")
        sys.exit(1)

    modify_fastq_headers(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

