import csv
from collections import defaultdict
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

# Initialize a dictionary to store the counts
counts = defaultdict(lambda: defaultdict(set))

# Read the CSV file and process it
with open(input_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        gRNA_name = row['gRNA_name'].strip()
        cell = row['cell'].strip()
        umi = row['umi'].strip()
        # Use a set to store unique UMIs
        counts[gRNA_name][cell].add(umi)

# Write the output to a TSV file
with open(output_file, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['gRNA_name', 'cell', 'count'])  # Write header
    for gRNA_name, cell_dict in counts.items():
        for cell, umi_set in cell_dict.items():
            writer.writerow([gRNA_name, cell, len(umi_set)])  # Write each row with the count of unique UMIs

