#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Summarize single-guide assignments into per-cell counts.")
    parser.add_argument("filtered", help="Filtered guide assignment table")
    parser.add_argument("output", help="Output TSV with aggregated counts")
    args = parser.parse_args()

    umi_sets = defaultdict(set)
    read_counts = defaultdict(int)

    with open(args.filtered, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell = row.get("Cell_Barcode")
            guide = row.get("Guide_Name") or row.get("SPACER1")
            umi = row.get("Molecular_Barcode")
            num_obs_raw = row.get("Num_Obs", "0")
            try:
                num_obs = int(num_obs_raw)
            except ValueError:
                try:
                    num_obs = int(float(num_obs_raw))
                except ValueError:
                    num_obs = 0
            key = (cell, guide)
            if umi:
                umi_sets[key].add(umi)
            read_counts[key] += num_obs

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["Cell_Barcode", "Guide_Name", "umi_count", "read_count"])
        for (cell, guide), reads in sorted(read_counts.items(), key=lambda kv: (kv[0][0] or "", kv[0][1] or "")):
            umi_count = len(umi_sets[(cell, guide)])
            writer.writerow([cell, guide, umi_count, reads])


if __name__ == "__main__":
    main()
