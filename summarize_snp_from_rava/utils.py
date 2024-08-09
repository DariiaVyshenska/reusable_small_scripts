import pandas as pd

def reorder_columns(summary_df):
  determ_col = [
    'SAMPLE_ID', 'NONSYNONYMOUS', 'SYNONYMOUS',  
    'COMPLEX', 'STOPGAIN', 'STOPLOSS', 'TOTAL_SNPS',
    ]
  unmatched_columns = [col for col in summary_df.columns if col not in determ_col]
  new_col_order = [col for col in determ_col if col in summary_df.columns] + unmatched_columns
  
  return summary_df[new_col_order]

def format_snp_table(snp_data):
  f_snp_data = snp_data[['Sample', 'Protein', 'AminoCorrect', 'NucleotideChange', 'AF', 'Syn']].copy()
  f_snp_data.columns = [
    'SAMPLE_ID', 'PROTEIN', 'AA_CHANGE', 'NT_CHANGE', 'ALLELE_FREQ', 'MUTATION_TYPE'
  ]
  f_snp_data['MUTATION_TYPE'] = f_snp_data['MUTATION_TYPE'].str.replace(' SNV', '')
  f_snp_data['SAMPLE_ID'] = f_snp_data['SAMPLE_ID'].str.replace('.fastq.gz', '')
  return f_snp_data