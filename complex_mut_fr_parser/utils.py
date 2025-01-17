def calc_freq(codon_counts, total_reads):
  return [(codon, round(count / total_reads, 4), count) for codon, count in codon_counts.items()]