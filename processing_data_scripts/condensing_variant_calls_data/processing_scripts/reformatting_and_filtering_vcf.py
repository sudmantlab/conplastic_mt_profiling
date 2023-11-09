#Input: all_muts.vcf
#Output: 
#1) somatic_mutations.vcf: all haplotype sites have been filtered from each sample and contains all heteroplasmic allele frequency information
#2) haplotype_mutations.vcf: this vcf contains only heteroplasmic frequency for the haplotype sites that define the mt haplotypes per strain

import pandas as pd
import numpy as np

#the input is the all_muts.vcf as specified in the rule
vcf = pd.read_csv(snakemake.input[0], sep = "\t")

vcf.reset_index(level=0, inplace=True)
vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ALLELE_INFO","FILE_NAME"]

#identifying each sample name from the filename
vcf["SAMPLE"] = vcf["FILE_NAME"].str.split(pat = "\.", expand = True)[0]

#creating identifiers for strain, tissue type, and age
vcf["STRAIN"] = vcf["SAMPLE"].str.split(pat = "_", expand = True)[0]
vcf["TISSUE"] = vcf["SAMPLE"].str.split(pat = "_", expand = True)[2]
vcf["AGE_BIN"] = vcf["SAMPLE"].str.split(pat = "_", expand = True)[1].str[0]
vcf["AGE_BIN"] = vcf["AGE_BIN"].map({"O":"OLD", "Y": "YOUNG", "M":"MID"})

#resetting coordinates of variants
vcf["START"] = vcf["POS"] - 1
vcf["END"] = vcf["POS"]

#we don't care for the filename column now
vcf_subset = vcf[["SAMPLE", "STRAIN", "TISSUE", "AGE_BIN", "START", "END", "REF", "ALT", "ALLELE_INFO"]]

#we want to add variant type information
def var_type(row):
    if len(row["REF"]) == 1 and len(row["ALT"]) == 1:
        type = "SNV"
    elif len(row["REF"]) > len(row["ALT"]):
        type = "DEL"
    else: 
        type = "INS"
    
    return type 

vcf_subset["VARIANT_TYPE"] = vcf_subset.apply(var_type, axis=1) 

#we parse the allele_info column now to separate the ref and alt across duplex molecules -- heteroplasmic frequencies
allele_info_split = vcf_subset["ALLELE_INFO"].str.split(pat=":", expand = True)

#REF,ALT
allele_info_split.columns = ["ALLELE_DEPTH", "READ_DEPTH_AT_POS", "FREQ", "COUNT_OF_NS"]
read_depth_count_of_ns = allele_info_split[["READ_DEPTH_AT_POS", "COUNT_OF_NS"]]

allele_depth_split = allele_info_split["ALLELE_DEPTH"].str.split(pat = ",", expand = True)
allele_depth_split.columns = ["REF_ALLELE_DEPTH", "ALT_ALLELE_DEPTH"]

freq_split = allele_info_split["FREQ"].str.split(pat = ",", expand = True)
freq_split.columns = ["REF_FREQ", "ALT_FREQ"]

vcf_allele_info = pd.concat([vcf_subset, allele_depth_split, freq_split, read_depth_count_of_ns], axis=1, ignore_index=True)

vcf_allele_info.columns = ["SAMPLE", "STRAIN", "TISSUE", "AGE_BIN", "START", "END", "REF", "ALT", "ALLELE_INFO",\
                          "VARIANT_TYPE", "REF_ALLELE_DEPTH", "ALT_ALLELE_DEPTH", "REF_FREQ", "ALT_FREQ","READ_DEPTH_AT_POS", "COUNT_OF_NS"]

#we don't care to keep the column with all heteroplasmic frequency info 
cleaned_vcf = vcf_allele_info[["SAMPLE", "STRAIN", "TISSUE", "AGE_BIN", "START", "END", "REF", "ALT",\
                          "VARIANT_TYPE","REF_ALLELE_DEPTH", "ALT_ALLELE_DEPTH", "REF_FREQ", "ALT_FREQ","READ_DEPTH_AT_POS", "COUNT_OF_NS"]]

#now we go into filtering out haplotype sites based on the conplastic strains
ALR = cleaned_vcf[cleaned_vcf["STRAIN"] == "ALR"]

POS = [4738, 9347, 9460]
REF = ["C", "G", "T"] 
ALT = ["A", "A", "C"]

ALR_hap_vars = pd.DataFrame({"START": POS, "REF": REF, "ALT": ALT})

keys = list(ALR_hap_vars.columns.values)
i_allvars = ALR.set_index(keys).index
i_hapvars = ALR_hap_vars.set_index(keys).index
ALR_haplotypes = ALR[i_allvars.isin(i_hapvars)]

FVB = cleaned_vcf[cleaned_vcf["STRAIN"] == "F"]

POS = [7777, 9460]
REF = ["G", "T"] 
ALT = ["T", "C"]

FVB_hap_vars = pd.DataFrame({"START": POS, "REF": REF, "ALT": ALT})

keys = list(FVB_hap_vars.columns.values)
i_allvars = FVB.set_index(keys).index
i_hapvars = FVB_hap_vars.set_index(keys).index
FVB_haplotypes = FVB[i_allvars.isin(i_hapvars)]

AKR = cleaned_vcf[cleaned_vcf["STRAIN"] == "AKR"]

POS = [9460]
REF = ["T"] 
ALT = ["C"]

AKR_hap_vars = pd.DataFrame({"START": POS, "REF": REF, "ALT": ALT})

keys = list(AKR_hap_vars.columns.values)
i_allvars = AKR.set_index(keys).index
i_hapvars = AKR_hap_vars.set_index(keys).index
AKR_haplotypes = AKR[i_allvars.isin(i_hapvars)]

NZB = cleaned_vcf[cleaned_vcf["STRAIN"] == "NZB"]

POS = [54, 1352, 1518, 1589, 1821, 2200, 2339, 2524, 2765, 2766, 2797, 2813, 2839, 2933, 3193, 3259, 3421, \
      3466, 3598, 3691, 3931, 4122, 4275, 4323, 4407, 4705, 4731, 4770, 4884, 4902, 5203, 5462, 5551, 5929, 6040, \
      6406, 6469, 6574, 6619, 6784, 7410, 7545, 7869, 8438, 8466, 8567, 8857, 8863, 9136, 9151, 9390, 9460, 9529, 9580, \
      9598, 9984, 10546, 10582, 10951, 11842, 11845, 11932, 12352, 12574, 12694, 12834, 12889, 13003, 13443,\
      13611, 13688, 13780, 13781, 13836, 13982, 14185, 14210, 14362, 14641, 14737, 15498, 15548, 15577, 15587, 15602,\
       15656, 15916, 16016, 16267, 16271]
REF = ["G","A","G","G","T","T","G","C","A","T","C","T","C","C","T",\
       "A","T","T","T","A","G","C","G","T","G","A","C","T","A","T","A","G",\
       "T","G","T","C","A","C","G","G","A","A","G","A","T","C", "T","C","A","T","A","T","C",\
       "C","A","G","C","A","C","G","C","A","C","T","A","T","A","G","C","T","C", \
       "A","T","A","A","T","G","A","G","C","T","C","A","C","C","T","C","A","A","T"] 

ALT = ["A","G","A","A","C","C","A","T","G","C","T","C","T","T","C","G","C","C","C",\
       "G","A","T","A","C","A","G","T","C","C","G","AG","A","C","A","C","T","G","T","A","A","G","G","A","G","C","T",\
      "C","T","G","C","G","C","T","T","G","A","T","G","A","A","T","C","T","A","G","C","G","A","T","C","T",\
      "G","C","G","G","C","A","G","A","T","A","T","T","T","T","C","T","C","G","C"]

NZB_hap_vars = pd.DataFrame({"START": POS, "REF": REF, "ALT": ALT})

keys = list(NZB_hap_vars.columns.values)
i_allvars = NZB.set_index(keys).index
i_hapvars = NZB_hap_vars.set_index(keys).index
NZB_haplotypes = NZB[i_allvars.isin(i_hapvars)]

#combining all entries which correspond to haplotype variants in the vcf file 
all_haplotypes = pd.concat([ALR_haplotypes, AKR_haplotypes, FVB_haplotypes, NZB_haplotypes], ignore_index=True, sort=False)

#filtering out all haplotype sites from the reformatted all_muts.vcf
keys = list(all_haplotypes.columns.values)
i_allvars = cleaned_vcf.set_index(keys).index
i_hapvars = all_haplotypes.set_index(keys).index

#splitting the somatic and haplotype vcfs
somatic_mutations_vcf = cleaned_vcf[~i_allvars.isin(i_hapvars)]
haplotype_mutations_vcf = cleaned_vcf[i_allvars.isin(i_hapvars)]

somatic_output_vcf = snakemake.output[0]
haplotype_output_vcf = snakemake.output[1]

#saving the files
somatic_mutations_vcf.to_csv(somatic_output_vcf, header=True, index=False, sep = "\t")
haplotype_mutations_vcf.to_csv(haplotype_output_vcf, header=True, index=False, sep = "\t")

