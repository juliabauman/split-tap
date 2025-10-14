import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def parse_args():
  parser = argparse.ArgumentParser(description="Filter guide assignments based on transcript-per-transcript thresholds.")
  parser.add_argument("input", help="Input TSV with UMI-level observations")
  parser.add_argument("output", help="Output TSV after TPT filtering")
  parser.add_argument("threshold", type=float, help="Minimum TPT to retain a transcript")
  parser.add_argument("histogram", nargs="?", help="Optional path for the TPT histogram TSV")
  return parser.parse_args()


args = parse_args()
infile = args.input
out_file = args.output
thresh = args.threshold

def calculate_tpt(umi_obs):

  # calculate the total number of reads per cell-umi barcode combination
  reads = umi_obs.groupby(['Cell_Barcode', 'Molecular_Barcode'])['Num_Obs'].sum()

  # convert to DataFrame with barcodes as columns and rename Num_obs column
  reads = pd.DataFrame(reads).reset_index()
  reads = reads.rename(columns = {'Num_Obs': 'total_reads'})

  # add total reads to umi_counts
  umi_obs = pd.merge(umi_obs, reads, how = 'left', on = ['Cell_Barcode', 'Molecular_Barcode'])

  # calculate the fraction of total reads for each gene-cell-umi combination, where total reads is
  # the sum of reads for that cell-umi combination. this metric is also called transcripts per transcript (tpt).
  umi_obs['tpt'] = np.divide(umi_obs['Num_Obs'], umi_obs['total_reads'])
  umi_tpt = umi_obs.drop(columns = 'total_reads')

  return(umi_tpt)

# compute tpt histogram from calculate_tpt() output
def compute_tpt_histogram(umi_tpt):
  
  # calculate number of bc-umi-gene combinations per tpt value
  tpt_hist = pd.value_counts(umi_tpt.tpt).to_frame().reset_index()
  
  # rename columns and sort according to tpt value
  tpt_hist.columns = ['tpt','transcripts']
  tpt_hist = tpt_hist.sort_values('tpt', ascending = False)
  
  return(tpt_hist)

# filter umi observations based on minimum tpt value
def filter_tpt(umi_tpt, tpt_threshold):

  # retain observations with tpt >= tpt_threshold
  filter_logical = umi_tpt['tpt'] >= tpt_threshold
  umi_filt = umi_tpt[filter_logical]

  print('  Removed ' + str(np.round(100 * (1 - np.mean(filter_logical)), 3)) + '% of transcripts')

  return(umi_filt)



# read umi observations input file
print('Reading input file...')
umi_obs = pd.read_csv(infile, sep = '\t')

print('Calculating TPT...')
umi_tpt = calculate_tpt(umi_obs)

# compute tpt histogram
tpt_hist = compute_tpt_histogram(umi_tpt)

# write histogram to text file
if args.histogram:
  hist_path = Path(args.histogram)
else:
  hist_path = Path(out_file).with_name("tpt_histogram.tsv")
hist_path.parent.mkdir(parents=True, exist_ok=True)
tpt_hist.to_csv(hist_path, sep = '\t', index = False)

# filter transcripts based on minimum tpt
print('Filtering for chimeric reads (min TPT = ' + str(thresh) + ')...')
umi_filt = filter_tpt(umi_tpt, tpt_threshold = thresh)
umi_filt.to_csv(out_file, sep = '\t', index = False)

print('Done!')
