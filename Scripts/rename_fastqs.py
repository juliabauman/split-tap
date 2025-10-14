import pandas as pd
import os
import shutil
from glob import glob


##run me from folder "admera"!!


# Load spreadsheet
df = pd.read_csv("../25095-01-QC-Summary-for-Seq-only.csv", header=0)

for _, row in df.iterrows():
    sample = row['Sample']                   # Provider's sample name
    customer_id = str(row['Customer_ID']).strip()  # Your sample name

    # Make a folder for each of your samples (e.g., A1_1)
    os.makedirs(customer_id, exist_ok=True)

    # Match the original FASTQ filenames from the provider
    pattern = f"{sample}_S*_L*_R*_001.fastq.gz"
    fastq_files = glob(pattern)

    for f in fastq_files:
        if "_R1_" in f:
            suffix = "R1"
        elif "_R2_" in f:
            suffix = "R2"
        else:
            continue  # Skip unrecognized file

        new_name = f"{customer_id}_{suffix}.fastq.gz"
        dest_path = os.path.join(customer_id, new_name)

        #print(f"[DRY RUN] Would move: {f} → {dest_path}")
        shutil.move(f, dest_path)
        print(f"Moved {f} → {dest_path}")
