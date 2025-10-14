import pandas as pd

df = pd.read_csv("Data/CROP_ids_dedup.txt", sep="\t")

# Group by Cell_Barcode and tRNA, then count UMIs
umi_counts = df.groupby(["Cell_Barcode", "tRNA"])["UMI"].nunique().reset_index(name="UMI_count")

# For each cell, find the tRNA with the most UMIs
dominant_tRNA = umi_counts.sort_values("UMI_count", ascending=False).drop_duplicates("Cell_Barcode")

# Save to file
dominant_tRNA.to_csv("Data/cell_trna_assignments.txt", sep="\t", index=False)

