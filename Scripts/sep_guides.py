import gzip
import sys

# Define the target sequence (must be uppercase)
TARGET_SEQ = "GACCCAGAAAAAAAGCACCGAC"

# Function to filter reads based on exact match in R1
def filter_reads(r1_file, r2_file, output_r1, output_r2):
    with gzip.open(r1_file, "rt") as r1, gzip.open(r2_file, "rt") as r2, \
         gzip.open(output_r1, "wt") as out_r1, gzip.open(output_r2, "wt") as out_r2:
        
        while True:
            # Read 4 lines per read (FASTQ format)
            r1_lines = [r1.readline().strip() for _ in range(4)]
            r2_lines = [r2.readline().strip() for _ in range(4)]
            
            # Stop when the end of the file is reached
            if not r1_lines[0]:
                break
            
            # Check if TARGET_SEQ is in the R1 sequence
            if TARGET_SEQ in r1_lines[1]:
                # Write the matching read pairs to output
                out_r1.write("\n".join(r1_lines) + "\n")
                out_r2.write("\n".join(r2_lines) + "\n")

r1_file=sys.argv[1]
r2_file=sys.argv[2]
output_r1=sys.argv[3]
output_r2=sys.argv[4]

filter_reads(r1_file, r2_file, output_r1, output_r2)
