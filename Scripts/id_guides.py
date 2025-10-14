import pandas as pd
import sys

def merge_spacer_files(txt_file_path, csv_file_path, output_file_path, ct_file):
    # Load the tab-separated TXT file
    df_txt = pd.read_csv(txt_file_path, sep="\t")

    # Load the Modified Spacer Pairs CSV file
    df_csv = pd.read_csv(csv_file_path)

    # Create the new column in the TXT file (spacer1 + spacer2)
    df_txt["spacer1_spacer2"] = df_txt["SPACER1"] + df_txt["SPACER2"]

    # Merge with the CSV file on the matching column
    df_merged = df_txt.merge(df_csv[["Seq1-Seq2_rc", "Guide1-Guide2_Name"]], 
                             left_on="spacer1_spacer2", 
                             right_on="Seq1-Seq2_rc", 
                             how="left")

    n_total = len(df_txt)
    n_merged = df_merged["Guide1-Guide2_Name"].notna().sum()
    n_unmatched = n_total - n_merged

    # Drop the extra Seq1_Seq2_rc column from the CSV merge
    df_merged.drop(columns=["Seq1-Seq2_rc","SPACER1","SPACER2"], inplace=True)

    df_dedup = df_merged.drop_duplicates()

    # Save the modified file
    df_dedup.to_csv(output_file_path, sep="\t", index=False)

    df_counts = df_dedup.groupby(['Cell_Barcode', 'Guide1-Guide2_Name'])['UMI'].count().reset_index()
    df_counts.rename(columns={'UMI': 'UMI_Count'}, inplace=True)
    df_counts.to_csv(ct_file, sep="\t", index=False)

    print(f"Total rows in input file: {n_total}")
    print(f"Matched to guide pair: {n_merged}")
    print(f"Unmatched rows (lost guide assignment): {n_unmatched} ({100 * n_unmatched / n_total:.2f}%)")

    print(f"File saved as {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python id_guides.py <input_txt_file> <csv_file> <output_txt_file> <counts_file>")
        sys.exit(1)

    input_txt = sys.argv[1]
    csv_file = sys.argv[2]
    output_txt = sys.argv[3]
    counts_file = sys.argv[4]

    merge_spacer_files(input_txt, csv_file, output_txt, counts_file)
