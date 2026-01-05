import pandas as pd
import sys
import numpy as np

infile = sys.argv[1]
out_file = sys.argv[2]
thresh = float(sys.argv[3])
out_hist = sys.argv[4]


#Calculate transcript per transcript (TPT) (Dixit A., 2016)
def calculate_tpt(umi_obs):
  #First: Get read counts of each guide-cell-UMI pairing
  for c in ['cell', 'umi', 'gRNA_name']:
        if c not in umi_obs.columns:
            raise ValueError(f"Missing required column: {c}")

  counts = (umi_obs.groupby(['cell','umi','gRNA_name'], as_index=False).size()
                    .rename(columns={'size':'Num_Obs'}))

  # 2. Compute total reads per (cell, umi)
  totals = (counts.groupby(["cell", "umi"])["Num_Obs"]
              .sum()
              .rename("total_reads")
              .reset_index())

  # 3. Merge totals back and compute TPT
  out = counts.merge(totals, on=["cell", "umi"], how="left")
  out["tpt"] = out["Num_Obs"] / out["total_reads"]

  return(out)


# Create tpt histogram from calculate_tpt() output
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
tpt_hist.to_csv(out_hist, sep = '\t', index = False)

# filter transcripts based on minimum tpt
print('Filtering for chimeric reads (min TPT = ' + str(thresh) + ')...')
umi_filt = filter_tpt(umi_tpt, tpt_threshold = thresh)
umi_filt.to_csv(out_file, sep = '\t', index = False)

print('Done!')
