import sys
import csv
import Bio.Align
from Bio import SeqIO
import gzip

#note: Use python3.9 versions of the modules

input_fastq = sys.argv[1]
output_txt = sys.argv[2]
guide_library_csv = sys.argv[3]
summary_data = []

print("Input FASTQ file: {0}".format(input_fastq))
print("Output TXT file: {0}".format(output_txt))
print("Guide library CSV: {0}".format(guide_library_csv))

# Load guide sequences from a CSV file into a dictionary
def load_guide_sequences(filename):
    guide_dict = {}
    with open(filename, mode='r') as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip header
        for rows in reader:
            guide_seq = rows[8]  # Assuming guide sequence is in the second column
            guide_name = rows[4]  # Assuming guide name is in the first column
            guide_dict[guide_seq] = guide_name
        keys = list(guide_dict.keys())[:5]
        for key in keys:
            print("{0}: {1}".format(key, guide_dict[key]))
    return guide_dict

guide_dict = load_guide_sequences(guide_library_csv)

start_seq = "GAAAGGACGAAACACC"
end_seq = "GTTTAAGAGCTATGCT"
max_mismatches = 1

def find_with_exact_matching(seq, pattern):
    return seq.find(pattern)  # Returns the position of the first match, or -1 if not found

def find_guide_name(spacer, guide_dict):
    return guide_dict.get(spacer, "No match found")

with open(output_txt, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['spacer', 'cell', 'umi', 'gRNA_name'])

    # Read the gzipped FASTQ file using gzip and SeqIO
    with gzip.open(input_fastq, "rt") as fastq_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            read_name = record.id
            read_seq = str(record.seq)

            cell_umi_parts = read_name.split("__")
            if len(cell_umi_parts) >= 5:
                cell = cell_umi_parts[3]
                umi = cell_umi_parts[4]
            else:
                continue  # Skip reads that don't have the expected format

            start_pos = find_with_exact_matching(read_seq, start_seq)
            end_pos = find_with_exact_matching(read_seq, end_seq)

            if start_pos != -1 and end_pos != -1:
                spacer_start = start_pos + len(start_seq)
                spacer = read_seq[spacer_start:end_pos] if spacer_start < end_pos and end_pos - spacer_start < 23 else ""
                gRNA_name = find_guide_name(spacer, guide_dict) if spacer else "No spacer found"
                # Write each row immediately
                writer.writerow([spacer, cell, umi, gRNA_name])
