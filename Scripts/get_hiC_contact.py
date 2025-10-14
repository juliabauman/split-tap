#Get Hi-C interaction frequencies for all E-G pairs on the same chromosome

##Run like python3 [scriptname.py] enhancer_pair_locations_for_HiC.csv hic_contacts_file desired_resolution
#maybe start with 5000 for resolution (what Andreas suggested)

import hicstraw
import pandas as pd
import sys

input_locations=sys.argv[1]
hic_file=sys.argv[2]
res=sys.argv[3]

# Define functions ---------------------------------------------------------------------------------

# get contact frequency between two bins
def get_int_freq_pair(bin1, bin2, matrix_chr):
    int_freq = matrix_chr.getRecords(int(bin1), int(bin1), int(bin2), int(bin2))
    if len(int_freq) > 0:
        return int_freq[0].counts
    else:
        return float("nan")

# general version for same-chromosome locus pairs
def get_pair_int_freqs_chr(etp, hic, chrom, hic_res, bin1_col, bin2_col, chr1_col, chr2_col):
    matrix_chr = hic.getMatrixZoomData(chrom, chrom, "observed", "NONE", "BP", hic_res)
    etp_chr = etp.loc[(etp[chr1_col] == chrom) & (etp[chr2_col] == chrom)].copy()
    etp_chr['hic_int_freq'] = etp_chr.apply(
        lambda x: get_int_freq_pair(x[bin1_col], x[bin2_col], matrix_chr),
        axis=1
    )
    return etp_chr

def get_pair_int_freqs(
    etp, hic, hic_res,
    coord1_start, coord1_end, chr1_col,
    coord2_start, coord2_end, chr2_col,
    bin1_col='pert1_start', bin2_col='pert2_start'
):
    # Compute center + bin for each locus
    etp['pert1_center'] = (etp[coord1_start] + etp[coord1_end]) / 2
    etp[bin1_col] = etp['pert1_center'] // hic_res * hic_res

    etp['pert2_center'] = (etp[coord2_start] + etp[coord2_end]) / 2
    etp[bin2_col] = etp['pert2_center'] // hic_res * hic_res

    # Compute frequencies for cis pairs
    print("Getting interaction frequencies for:")
    etp_hic = []
    for chrom in etp[chr1_col].unique():
        print(f"  {chrom}")
        etp_hic.append(
            get_pair_int_freqs_chr(etp, hic, chrom, hic_res, bin1_col, bin2_col, chr1_col, chr2_col)
        )
    etp_hic = pd.concat(etp_hic)

    # Fill in trans-pairs with NaN
    etp_trans = etp.loc[etp[chr1_col] != etp[chr2_col]].copy()
    etp_trans['hic_int_freq'] = float("nan")
    etp_hic = pd.concat([etp_hic, etp_trans])

    return etp_hic


### Get contact frequencies for E-G pairs ------------------------------------------------------------

# load ETP table
etp = pd.read_csv(input_locations)

# open connection to .hic file
hic = hicstraw.HiCFile(hic_file)

# get interaction frequencies for all E-G pairs in ETP table and write output to file
etp = get_pair_int_freqs(
    etp, hic, hic_res=int(res),
    coord1_start='pert1_start', coord1_end='pert1_end', chr1_col='pert1_chr',
    coord2_start='pert2_start', coord2_end='pert2_end', chr2_col='pert2_chr'
)

etp = etp.drop(columns=['pert1_center', 'pert1_start', 'pert2_end', 'pert2_start', 'pert1_end', 'pert2_center'])
etp.to_csv("hic_contact_frequencies.csv", index = False)

