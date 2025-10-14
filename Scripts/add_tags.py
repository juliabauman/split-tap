
import pysam
import sys

input_bam = sys.argv[1]
output_bam = sys.argv[2]

with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
    for read in infile:
        # Extract parts of the read name
        parts = read.query_name.split("__")
        if len(parts) > 3:
            cell_barcode = parts[3]
            umi = parts[4]
            # Add the CB:Z and UB:Z tags
            read.set_tag("CB", cell_barcode, value_type="Z")
            read.set_tag("UB", umi, value_type="Z")
        # Write the modified read to the output file
        outfile.write(read)

