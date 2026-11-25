#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict, Counter


def merge_spacer_files(input_txt, csv_file, output_txt, counts_file, nomatch_file):
    # -----------------------------
    # Load spacer-pair → guide map
    # -----------------------------
    spacer_to_guide = {}
    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        required_cols = {"Seq1-Seq2_rc", "Guide1-Guide2_Name"}
        if not required_cols.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"CSV missing required columns. Expected {required_cols}, "
                f"found {set(reader.fieldnames or [])}"
            )
        for row in reader:
            spacer_to_guide[row["Seq1-Seq2_rc"].strip()] = row["Guide1-Guide2_Name"].strip()

    print(f"Loaded {len(spacer_to_guide):,} spacer-pair mappings", flush=True)

    # -----------------------------
    # Streaming processing
    # -----------------------------
    counts = defaultdict(int)          # (cell, guide) -> unique UMI count
    seen_umis = defaultdict(set)       # (cell, guide) -> set of UMIs seen
    nomatch = Counter()                # spacer1_spacer2 -> count

    n_total = 0
    n_matched_rows = 0
    n_unmatched_rows = 0
    n_unique_written = 0
    n_bad_rows = 0

    with open(input_txt, newline="") as fin, open(output_txt, "w", newline="") as fout:
        tsv_reader = csv.reader(fin, delimiter="\t")
        header = next(tsv_reader, None)
        if header is None:
            raise ValueError("Input TXT is empty (no header row).")

        # Make sure tRNA exists
        if "tRNA" not in header:
            raise ValueError(f"Input TXT missing required column: tRNA. Columns: {header}")

        # Build indices (safe because we now trust csv parsing)
        idx = {c: i for i, c in enumerate(header)}

        required = ["cell", "umi", "SPACER1", "SPACER2", "tRNA"]
        missing = [c for c in required if c not in idx]
        if missing:
            raise ValueError(f"Input TXT missing required columns: {missing}")

        # Output header expected downstream
        fout.write("cell\tumi\ttRNA\tGuide1-Guide2_Name\n")

        for i, fields in enumerate(tsv_reader, start=2):  # start=2 because header is line 1
            n_total += 1

            # If row has wrong number of columns, skip + count it (prevents index drift)
            if len(fields) != len(header):
                n_bad_rows += 1
                # If you prefer to hard-fail instead, replace with: raise ValueError(...)
                continue

            spacer1 = fields[idx["SPACER1"]].strip()
            spacer2 = fields[idx["SPACER2"]].strip()
            cell    = fields[idx["cell"]].strip()
            umi     = fields[idx["umi"]].strip()
            trna    = fields[idx["tRNA"]].strip()

            sp12 = spacer1 + spacer2  # same as your pandas script
            guide = spacer_to_guide.get(sp12)

            if guide is None:
                n_unmatched_rows += 1
                nomatch[sp12] += 1
                continue

            n_matched_rows += 1

            # Unique UMI per (cell, guide)
            cg = (cell, guide)
            if umi in seen_umis[cg]:
                continue
            seen_umis[cg].add(umi)

            counts[cg] += 1
            fout.write(f"{cell}\t{umi}\t{trna}\t{guide}\n")
            n_unique_written += 1

            if n_total % 5_000_000 == 0:
                print(
                    f"Processed {n_total:,} data rows | matched: {n_matched_rows:,} | "
                    f"unmatched: {n_unmatched_rows:,} | unique written: {n_unique_written:,} | "
                    f"bad_rows_skipped: {n_bad_rows:,}",
                    flush=True
                )

    # -----------------------------
    # Write counts file
    # -----------------------------
    with open(counts_file, "w", newline="") as f:
        f.write("cell\tGuide1-Guide2_Name\tUMI_Count\n")
        for (cell, guide), n in counts.items():
            f.write(f"{cell}\t{guide}\t{n}\n")

    # -----------------------------
    # Write nonmatch summary
    # -----------------------------
    with open(nomatch_file, "w", newline="") as f:
        f.write("spacer1_spacer2\tcount\n")
        for sp12, c in nomatch.most_common():
            f.write(f"{sp12}\t{c}\n")

    print(f"Total input data rows read: {n_total:,}", flush=True)
    print(f"Bad/malformed rows skipped (wrong # columns): {n_bad_rows:,}", flush=True)
    print(f"Matched rows (had guide mapping): {n_matched_rows:,}", flush=True)
    print(f"Unmatched rows: {n_unmatched_rows:,} ({100*n_unmatched_rows/max(n_total,1):.2f}%)", flush=True)
    print(f"Unique UMI rows written to output: {n_unique_written:,}", flush=True)
    print(f"Wrote: {output_txt}", flush=True)
    print(f"Wrote: {counts_file}", flush=True)
    print(f"Wrote: {nomatch_file}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python id_guides.py <input_txt_file> <csv_file> <output_txt_file> <counts_file> <nomatch_file>")
        sys.exit(1)

    input_txt   = sys.argv[1]
    csv_file    = sys.argv[2]
    output_txt  = sys.argv[3]
    counts_file = sys.argv[4]
    nomatch_file= sys.argv[5]

    merge_spacer_files(input_txt, csv_file, output_txt, counts_file, nomatch_file)
