#This file takes as input the all_mutfreq.csv file and outputs two files: 1) overall mutation frequency file, with a calculated mutfreq for each sample taking all muts into consideration and 2) a mutfreq file where the mutfreq is split per variant type.

import pandas as pd
import numpy as np

mut_freq = pd.read_csv(snakemake.input[0], sep = ",")
mut_freq.rename(columns = {"#SAMPLE": "SAMPLE"}, inplace = True)

#creating fields to denote strain, tissue, and age
mut_freq["STRAIN"] = mut_freq["SAMPLE"].str.split(pat = "_", expand = True)[0]
mut_freq["TISSUE"] = mut_freq["SAMPLE"].str.split(pat = "_", expand = True)[2]
mut_freq["AGE_BIN"] = mut_freq["SAMPLE"].str.split(pat = "_", expand = True)[1].str[0]
mut_freq["AGE_BIN"] = mut_freq["AGE_BIN"].map({"O":"OLD", "Y": "YOUNG", "M":"MID"})

mut_freq_by_type = mut_freq[["SAMPLE","MUTATION_TYPE", "MUTATION_CLASS", "COUNT", "DENOMINATOR","FREQUENCY", "STRAIN","TISSUE", "AGE_BIN"]]

mut_freq_by_type_output = snakemake.output[0]
mut_freq_by_type.to_csv(mut_freq_by_type_output, header=True, index=False, sep = "\t")

#summing the counts across all variant types 
#we need to filter for the entries that refer to aggregated mutation frequencies and sequencing denominator
total_mut_freq = mut_freq[mut_freq["MUTATION_TYPE"].str.contains("Total")]
overall_mut_freq = total_mut_freq.groupby(["SAMPLE", "DENOMINATOR", "STRAIN", "TISSUE", "AGE_BIN"]).COUNT.sum().reset_index()
overall_mut_freq["MUT_FREQ"] = overall_mut_freq["COUNT"]/overall_mut_freq["DENOMINATOR"]

overall_mut_freq_output = snakemake.output[1]
overall_mut_freq.to_csv(overall_mut_freq_output, header=True, index=False, sep = "\t")

