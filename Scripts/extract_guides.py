import gzip
import sys
import pandas as pd

def extract_info(fastq_file, output_txt, summary_txt):
    """Extracts Cell Barcode, UMI, SPACER1, and SPACER2 from FASTQ headers and saves to a text file."""
    total_reads = 0
    skipped_reads = 0
    malformed_headers = 0

    with gzip.open(fastq_file, "rt") as fq, open(output_txt, "w") as out_txt, open(summary_txt, "w") as summary_out:
        # Write header for output file
        out_txt.write("Cell_Barcode\tUMI\tSPACER1\tSPACER2\ttRNA\n")

        while True:
            # Read the FASTQ header (every 4th line is a header)
            header = fq.readline().strip()
            fq.readline()  # Skip sequence
            fq.readline()  # Skip "+"
            fq.readline()  # Skip quality score

            # Stop when EOF is reached
            if not header:
                break

            total_reads += 1  # Count total reads

            # Extract Cell Barcode (everything before UMI)
            parts = header.split("__")
            if len(parts) < 5:
                print(f"Skipping malformed header: {header}")
                skipped_reads += 1
                malformed_headers += 1
                continue
            
            cell_barcode = parts[3]
            umi = parts[-1].split(",")[0]  # UMI is before the SPACER1 tag

            # Extract SPACER1 and SPACER2
            spacer1 = ""
            spacer2 = "NOSEQ"
            if ",SPACER1=" in header:
                spacer1 = header.split(",SPACER1=")[-1].split(":")[0]
            if ":SPACER2=" in header:
                spacer2 = header.split(":SPACER2=")[-1].split(":")[0]
            if "tRNA=" in header:
                trna = header.split("tRNA=")[-1]

            # Skip writing the row if one or more spacer is NOSEQ
            if not spacer1 or spacer2 == "NOSEQ":
                skipped_reads += 1
                continue

            out_txt.write(f"{cell_barcode}\t{umi}\t{spacer1}\t{spacer2}\t{trna}\n")

        # Write summary
        summary_out.write(f"Total reads processed: {total_reads}\n")
        summary_out.write(f"Reads skipped (SPACER1 or SPACER2 were NOSEQ): {skipped_reads}\n")
        summary_out.write(f"Valid reads written to output: {total_reads - skipped_reads}\n")
        summary_out.write(f"Malformed headers skipped: {malformed_headers}\n")


    print(f"✅ Extracted info saved to: {output_txt}")
    print(f"📊 Summary saved to: {summary_txt}")


def add_num_obs_column(infile, outfile_with_counts):
    """
    Reads input file, counts occurrences of each unique (Cell_Barcode, UMI, SPACER1, SPACER2, tRNA) combo,
    and saves a new file with a Num_Obs column.
    """
    df = pd.read_csv(infile, sep="\t")
    df_counts = df.groupby(['Cell_Barcode', 'UMI', 'SPACER1', 'SPACER2', 'tRNA']).size().reset_index(name='Num_Obs')
    df_counts.rename(columns={'UMI': 'Molecular_Barcode'}, inplace=True)
    df_counts.to_csv(outfile_with_counts, sep="\t", index=False)
    print(f"File with Num_Obs saved as {outfile_with_counts}")
    return df_counts


# Command-line execution
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_spacer_info.py input.fastq.gz output.txt summary.txt read_count.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    summary_file = sys.argv[3]
    CROP_read_ct = sys.argv[4]

    extract_info(input_file, output_file, summary_file)
    add_num_obs_column(output_file,CROP_read_ct)
