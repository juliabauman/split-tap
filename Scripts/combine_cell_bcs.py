import gzip
import sys

def modify_fastq_headers(input_fastq, output_fastq, index):
    with gzip.open(input_fastq, 'rt') as infile, gzip.open(output_fastq, 'wt') as outfile:
        while True:
            # Read four lines per read (FASTQ format)
            header = infile.readline().strip()
            if not header:
                break  # Stop at EOF
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Parse the header
            parts = header.split('__')
            
            if len(parts) < 5:
                print(f"Skipping malformed header: {header}")
                continue

            cell_barcode = parts[3]  # Extract original cell barcode
            umi = parts[4]
            trash = parts[5]  # Throw away
            
            # Append i7 index to the cell barcode
            i7_index=index
            new_cell_barcode = f"{cell_barcode}_{i7_index}"

            # Reconstruct the new header
            new_parts = parts[:3] + [new_cell_barcode] + [umi]
            new_header = '__'.join(new_parts)

            # Write the modified read to the output file
            outfile.write(f"{new_header}\n{sequence}\n{plus}\n{quality}\n")

    print(f"✅ Modified headers written to: {output_fastq}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python modify_fastq_headers.py input.fastq.gz output.fastq.gz index_seq")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    index_seq = sys.argv[3]
    
    modify_fastq_headers(input_file, output_file, index_seq)
