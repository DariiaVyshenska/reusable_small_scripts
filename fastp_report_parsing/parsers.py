import json

def parse_fastp_stats(json_file_path):
  with open(json_file_path, 'r') as file:
    data = json.load(file)
  
  total_reads = data["summary"]["before_filtering"]["total_reads"]
  duplication_rate = data["duplication"]["rate"]
  insert_size_peak = data["insert_size"]["peak"]
  insert_size_unknown = data["insert_size"]["unknown"]
  gc_content_before = data["summary"]["before_filtering"]["gc_content"]
  gc_content_after = data["summary"]["after_filtering"]["gc_content"]

  passed_filter_reads = data["filtering_result"]["passed_filter_reads"]
  low_quality_reads = data["filtering_result"]["low_quality_reads"]
  too_many_N_reads = data["filtering_result"]["too_many_N_reads"]
  too_short_reads = data["filtering_result"]["too_short_reads"]


  passed_filter_ratio = (passed_filter_reads / total_reads) * 100
  low_quality_ratio = (low_quality_reads / total_reads) * 100
  too_many_N_ratio = (too_many_N_reads / total_reads) * 100
  too_short_ratio = (too_short_reads / total_reads) * 100
  
  r1_len = data["summary"]["after_filtering"]["read1_mean_length"]
  r2_len = data["summary"]["after_filtering"]["read2_mean_length"]
  r1_len_before = data["summary"]["before_filtering"]["read1_mean_length"]
  r2_len_before = data["summary"]["before_filtering"]["read2_mean_length"]
  
  read_len_issues = (r1_len < (r1_len_before / 2) or r2_len < (r2_len_before / 2))

  
  if (read_len_issues):
    print(f"ATTENTION! File {json_file_path} has read length issues")

  return [
    total_reads, r1_len, r2_len, duplication_rate,
    insert_size_peak, insert_size_unknown,
    gc_content_before, gc_content_after,
    passed_filter_reads, passed_filter_ratio,
    low_quality_reads, low_quality_ratio,
    too_many_N_reads, too_many_N_ratio,
    too_short_reads, too_short_ratio,
  ]
