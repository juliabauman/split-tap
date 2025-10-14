#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Summarize dual-guide assignments into per-cell counts.")
    parser.add_argument("filtered", help="Filtered guide assignment table")
    parser.add_argument("output", help="Output TSV with aggregated counts")
    args = parser.parse_args()

    umi_sets = defaultdict(set)
    read_counts = defaultdict(int)

    with open(args.filtered, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell = row.get("Cell_Barcode")
            spacer1 = row.get("SPACER1")
            spacer2 = row.get("SPACER2")
            trna = row.get("tRNA", "")
            umi = row.get("Molecular_Barcode")
            num_obs_raw = row.get("Num_Obs", "0")
            try:
                num_obs = int(num_obs_raw)
            except ValueError:
                try:
                    num_obs = int(float(num_obs_raw))
                except ValueError:
                    num_obs = 0
            key = (cell, spacer1, spacer2, trna)
            if umi:
                umi_sets[key].add(umi)
            read_counts[key] += num_obs

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow([
            "Cell_Barcode",
            "SPACER1",
            "SPACER2",
            "tRNA",
            "umi_count",
            "read_count",
        ])
        for (cell, spacer1, spacer2, trna), reads in sorted(read_counts.items(), key=lambda kv: tuple((kv[0][i] or "") for i in range(4))):
            umi_count = len(umi_sets[(cell, spacer1, spacer2, trna)])
            writer.writerow([cell, spacer1, spacer2, trna, umi_count, reads])


if __name__ == "__main__":
    main()
