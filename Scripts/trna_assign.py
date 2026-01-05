import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t")

# Group by Cell_Barcode and tRNA, then count UMIs
umi_counts = df.groupby(["cell", "tRNA"])["umi"].nunique().reset_index(name="UMI_count")

# For each cell, find the tRNA with the most UMIs
dominant_tRNA = umi_counts.sort_values("UMI_count", ascending=False).drop_duplicates("cell")

# Save to file
dominant_tRNA.to_csv(output_file, sep="\t", index=False)

