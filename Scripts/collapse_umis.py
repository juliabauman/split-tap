import pandas as pd
import sys
from collections import Counter

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def collapse_by_hamming(sequences, counts, max_dist=1):
    seq_counts = Counter()
    for seq, count in zip(sequences, counts):
        seq_counts[seq] += count

    collapsed = set()
    result = {}

    for seq in sorted(seq_counts, key=seq_counts.get, reverse=True):
        if seq in collapsed:
            continue
        group = [s for s in seq_counts if s not in collapsed and hamming(s, seq) <= max_dist]
        for g in group:
            result[g] = seq
            collapsed.add(g)
    return result

input_file = sys.argv[1]
output_file = sys.argv[2]
umi_corrections_file = sys.argv[3]
spacer_corrections_file = sys.argv[4]

df = pd.read_csv(input_file, sep='\t')

# ── Stage 1: Collapse UMIs ────────────────────────────────────────────────────
umi_corrections = []
umi_collapsed_records = []

group_cols = ['Cell_Barcode', 'SPACER1', 'SPACER2']
for group_key, group_df in df.groupby(group_cols):
    umis = group_df['Molecular_Barcode'].tolist()
    counts = group_df['Num_Obs'].tolist()
    collapse_map = collapse_by_hamming(umis, counts)

    n_before = len(set(umis))
    n_after = len(set(collapse_map.values()))
    umi_corrections.append({
        'Cell_Barcode': group_key[0],
        'SPACER1': group_key[1],
        'SPACER2': group_key[2],
        'n_before': n_before,
        'n_after': n_after,
        'n_corrected': n_before - n_after
    })

    group_df['collapsed_umi'] = group_df['Molecular_Barcode'].map(collapse_map)
    collapsed = group_df.groupby(['collapsed_umi'], as_index=False).agg({'Num_Obs': 'sum'})
    for col in group_cols:
        collapsed[col] = group_key[group_cols.index(col)]
    umi_collapsed_records.append(collapsed)

umi_df = pd.concat(umi_collapsed_records, ignore_index=True)
umi_df = umi_df.rename(columns={'collapsed_umi': 'Molecular_Barcode'})

# ── Stage 2: Collapse SPACER1/SPACER2 ─────────────────────────────────────────
spacer_corrections = []
spacer_collapsed_records = []

group_cols = ['Cell_Barcode', 'Molecular_Barcode']
for group_key, group_df in umi_df.groupby(group_cols):
    spacer1s = group_df['SPACER1'].tolist()
    counts1 = group_df['Num_Obs'].tolist()
    collapse_map1 = collapse_by_hamming(spacer1s, counts1)

    spacer2s = group_df['SPACER2'].tolist()
    counts2 = group_df['Num_Obs'].tolist()
    collapse_map2 = collapse_by_hamming(spacer2s, counts2)

    spacer_corrections.append({
        'Cell_Barcode': group_key[0],
        'Molecular_Barcode': group_key[1],
        'n_spacer1_corrected': len(set(spacer1s)) - len(set(collapse_map1.values())),
        'n_spacer2_corrected': len(set(spacer2s)) - len(set(collapse_map2.values()))
    })

    group_df['collapsed_spacer1'] = group_df['SPACER1'].map(collapse_map1)
    group_df['collapsed_spacer2'] = group_df['SPACER2'].map(collapse_map2)

    collapsed = group_df.groupby(['collapsed_spacer1', 'collapsed_spacer2'], as_index=False).agg({'Num_Obs': 'sum'})
    for col in group_cols:
        collapsed[col] = group_key[group_cols.index(col)]
    spacer_collapsed_records.append(collapsed)

final_df = pd.concat(spacer_collapsed_records, ignore_index=True)
final_df = final_df.rename(columns={'collapsed_spacer1': 'SPACER1', 'collapsed_spacer2': 'SPACER2'})
final_df.to_csv(output_file, sep='\t', index=False)

pd.DataFrame(umi_corrections).to_csv(umi_corrections_file, sep='\t', index=False)
pd.DataFrame(spacer_corrections).to_csv(spacer_corrections_file, sep='\t', index=False)

