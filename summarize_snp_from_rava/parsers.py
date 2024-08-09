import pandas as pd
from utils import reorder_columns

def get_haf_info(snp_data, af_threshold):
  haf_snp_info = snp_data[(snp_data['ALLELE_FREQ'] >= af_threshold)]
  return haf_snp_info

def get_snp_summary(snp_info):
  summary_df = (
      snp_info[['SAMPLE_ID', 'MUTATION_TYPE']]
      .groupby(['SAMPLE_ID', 'MUTATION_TYPE'])
      .size()
      .unstack(fill_value=0)
      .reset_index()
    )
  summary_df['TOTAL_SNPS'] = summary_df.iloc[:, 1:].sum(axis=1)
  summary_df.columns = summary_df.columns.str.upper()
  return reorder_columns(summary_df)