#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO


def load_guide_reference(path):
    guide_map = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            try:
                guide_name = row[4].strip()
                guide_seq = row[8].strip()
            except IndexError:
                raise ValueError("Guide reference must have name in column 5 and sequence in column 9")
            if guide_seq:
                guide_map[guide_seq] = guide_name
    return guide_map


def extract_spacer(sequence, start_seq, end_seq):
    start_idx = sequence.find(start_seq)
    if start_idx == -1:
        return ""
    start_idx += len(start_seq)
    end_idx = sequence.find(end_seq, start_idx)
    if end_idx == -1 or end_idx <= start_idx:
        return ""
    return sequence[start_idx:end_idx]


def main():
    parser = argparse.ArgumentParser(description="Process single-guide CROP-seq reads into aggregated guide assignments.")
    parser.add_argument("fastq", help="Input FASTQ file with cell barcode/UMI annotated headers")
    parser.add_argument("reference", help="CSV file describing guide sequences (guide name column 5, sequence column 9)")
    parser.add_argument("assignments", help="Output TSV with per-read assignments")
    parser.add_argument("aggregated", help="Output TSV with aggregated per cell/UMI counts")
    parser.add_argument("--start-seq", default="GAAAGGACGAAACACC", help="Sequence immediately upstream of the guide spacer")
    parser.add_argument("--end-seq", default="GTTTAAGAGCTATGCT", help="Sequence immediately downstream of the guide spacer")

    args = parser.parse_args()

    guide_map = load_guide_reference(args.reference)
    assignments_path = Path(args.assignments)
    assignments_path.parent.mkdir(parents=True, exist_ok=True)
    aggregated_path = Path(args.aggregated)
    aggregated_path.parent.mkdir(parents=True, exist_ok=True)

    assignments = []
    counts = Counter()

    with open(args.fastq, "rt") as fastq_handle:
        for record in SeqIO.parse(fastq_handle, "fastq"):
            read_name = record.id
            parts = read_name.split("__")
            if len(parts) < 5:
                continue
            cell = parts[3]
            umi = parts[4]
            spacer = extract_spacer(str(record.seq), args.start_seq, args.end_seq)
            if not spacer:
                continue
            guide_name = guide_map.get(spacer, "No match found")
            assignments.append((spacer, cell, umi, guide_name))
            counts[(cell, umi, spacer, guide_name)] += 1

    with assignments_path.open("w", newline="") as assign_handle:
        writer = csv.writer(assign_handle, delimiter="\t")
        writer.writerow(["spacer", "cell", "umi", "guide_name"])
        writer.writerows(assignments)

    with aggregated_path.open("w", newline="") as agg_handle:
        writer = csv.writer(agg_handle, delimiter="\t")
        writer.writerow([
            "Cell_Barcode",
            "Molecular_Barcode",
            "SPACER1",
            "SPACER2",
            "Guide_Name",
            "Num_Obs",
        ])
        for (cell, umi, spacer, guide_name), num_obs in counts.items():
            writer.writerow([cell, umi, spacer, "NA", guide_name, num_obs])


if __name__ == "__main__":
    main()
