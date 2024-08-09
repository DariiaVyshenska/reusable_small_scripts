from parsers import get_snp_summary

def dump_all_snp_info(f_snp_data, output_path):
  all_snp_summary = get_snp_summary(f_snp_data)
  output_file_path_snp_summary = f"{output_path}/snp_summary.csv"
  all_snp_summary.to_csv(output_file_path_snp_summary, index=False)

def dump_haf_info(haf_snp_info, af_threshold, output_path):
  output_file_path_snp_list = f"{output_path}/af_greater{af_threshold}_snp_info.csv"
  haf_snp_info.to_csv(output_file_path_snp_list, index=False)
  
def dump_haf_summary(haf_snp_info, af_threshold, output_path):
  haf_snp_summary = get_snp_summary(haf_snp_info)
  output_file_path_haf_snp_summary = f"{output_path}/af_greater{af_threshold}_snp_summary.csv"
  haf_snp_summary.to_csv(output_file_path_haf_snp_summary, index=False)