import pandas as pd 
import numpy as np
import io
import os
import math
from Bio import SeqIO

#import the mm10 reference sequence
ref_seq = SeqIO.read(snakemake.input[0], "fasta")
ref_seq = ref_seq.upper()

#create a dictionary of the SNV haplotypes -- will allow us to create our reference genomes for all strains
#file needed to make the haplotype dictionary 
haplotype_info_file = snakemake.input[1]
haplotype_info = pd.read_csv(haplotype_info_file,  sep = "\t")

#filter out indels from the haplotype info (won't be included in our dictionary to mutate). Note: these indels are in noncoding regions
haplotype_info = haplotype_info[haplotype_info["VARIANT_TYPE"] == "SNV"]

#keep only one record of the haplotypes per strain 
haplotype_info = haplotype_info.drop_duplicates(subset = ["STRAIN", "START"], keep = "first")[["STRAIN", \
                                                                                         "START", "REF", "ALT"]]
strains = list(haplotype_info["STRAIN"].unique())

#create a dictionary with the structure {strain: [START, REF, ALT]}
haplotype_dict = {}
for strain in strains:
    #for a given strain, take the position, ref, and alt alleles of all the haploptypes & place into one list of lists 
    haplotypes = haplotype_info[haplotype_info["STRAIN"] == strain]
    haplotypes = haplotypes[["START", "REF", "ALT"]].values.tolist()
    #add the list of lists as a value to the strain keys
    haplotype_dict[strain] = haplotypes

#These are the two main functions we will be using
#get_ref_genome -- inputs: strain name & haplotype_dict; output: a seq object that contains the reference genome for the given haplotype strain
#get_all_possible_muts --  inputs: strain & strain_ref_genome; output: a dataframe with the following columns STRAIN POS (IN 0-INDEX) REF ALT. The dataframe contains all possible mutations at a given position in the genome Note: we begin the combinations at ND1 and finish at CYTB while containing all non-protein coding regions in between these genes

def get_ref_genome(strain, haplotype_dict):
    #we need to convert from a seq to a mutable object in order to mutate our sequence
    strain_ref_genome = ref_seq.seq.tomutable()

    haplotypes = haplotype_dict[strain]
    
    #we loop through each haplotype and mutate the ref sequence from mm10 accordingly
    for index in range(0, len(haplotypes)):
        site = haplotypes[index][0]
        allele = haplotypes[index][2]
        strain_ref_genome[site] = allele
    
    #convert back to a seq object in order to access all the methods 
    strain_ref_genome = strain_ref_genome.toseq()
    
    return strain_ref_genome


def get_all_possible_muts(strain, strain_ref_genome):
    
    #convert the sequence object into a dataframe with columns POS and REF (0-indexed position)
    strain_ref_df = pd.DataFrame([index, letter] for index,letter in enumerate(strain_ref_genome))
    strain_ref_df.columns = ["POS", "REF"]

    #create a dataframe that has the 4-nucleotides with the column name ALT
    nucleotide = pd.DataFrame(['A','G','C','T'])
    nucleotide.columns = ["ALT"]

    #the neat trick: give both dataframes a key that we will merge on 
    nucleotide['key'] = 1
    strain_ref_df['key'] = 1

    #the merge allows us to form all possible combinations of mutations at each site
    all_combos_df = pd.merge(strain_ref_df, nucleotide, on = 'key').drop('key', axis=1)

    #delete rows such that the reference allele = the alternative allele
    all_possible_muts_df = all_combos_df[all_combos_df["REF"] != all_combos_df["ALT"]]
    #keep entries between ND1 % CYTB, taking out ~3,000 positions we don't need for annotation 
    all_possible_muts_df = all_possible_muts_df[(all_possible_muts_df["POS"] >= 2750) &\
                                                (all_possible_muts_df["POS"] <= 15288)].reset_index(drop = True)
    all_possible_muts_df.insert(0, "STRAIN", strain)
    
    return all_possible_muts_df

#We call our main functions in this chunk of code and create a dataframe that contains all possible mutations for all strains  into one large dataframe (dim. 188085 rows by 4 columns)
all_possible_muts_df = pd.DataFrame()

#iterate through all 5 strains to create the all possible muts dataframe and append the resulting dataframe 
#to our final all_possible_muts_df 
for strain in ["AKR", "ALR", "B6", "F", "NZB"]:
    if strain == "B6":
        strain_ref_genome = ref_seq
    else:
        strain_ref_genome = get_ref_genome(strain, haplotype_dict)

    all_possible_muts_df = all_possible_muts_df.append(get_all_possible_muts(strain, strain_ref_genome),\
                                                       ignore_index = True)

output_file = snakemake.output[0]
all_possible_muts_df.to_csv(output_file, sep = "\t", index = False)