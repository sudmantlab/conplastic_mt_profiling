#IMS 08-21-2023
#input:
#output:
#parameters:
#important dependencies:

import pandas as pd
import numpy as np
#the following imports are from biopython - gives us access to reading in a fasta, Seq, SeqRecord, and MutableSeq functions
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
#the following 4 packages are part of base python
import sys
import io
import os
import argparse

#adding user defined inputs for our simulator
parser = argparse.ArgumentParser()
#parameters defining our experimental unit 
parser.add_argument("-strain", help = "conplastic strain being processed, IN ALL CAPS")
parser.add_argument("-tissue", help = "tissue being processed, Sentence case (eg Brain)")
parser.add_argument("-age", help = "age bin for samples, IN ALL CAPS")
parser.add_argument("-bin_freq", help = "the bin frequency being processed")

parser.add_argument("-sim_num", help = "Number of sims you wish to generate")

#input files needed to run this program 
parser.add_argument("-fasta_file_path", help = "File path to the corresponding fasta file")
parser.add_argument("-mut_count_file_path", help = "File path to the mut count file")
parser.add_argument("-mut_prop_file_path", help = "File path to the mut prop file")
parser.add_argument("-haplotype_file_path", help = "File path to a file containing information about mouse haplotypes that distinguish it from mm10")
parser.add_argument("-annotated_variants_file_path", help = "File path to a file containing annotations for all possible variants")
args = parser.parse_args()

#assigning our arguments to variables
strain = args.strain
tissue = args.tissue
age = args.age
bin_freq = args.bin_freq

sim_num = int(args.sim_num)

#load our mut_prop and mut_total_count files
mut_prop_file = args.mut_prop_file_path 
mut_total_count_file = args.mut_count_file_path 

mut_type_prop = pd.read_csv(mut_prop_file, sep = "\t")
mut_total_count = pd.read_csv(mut_total_count_file, sep = "\t")

#some conditions x bin combinations have 0 counts -- we want to export an empty file in that case (to keep as record) 
check_flag = mut_total_count[(mut_total_count["STRAIN"] == strain) &\
                                  (mut_total_count["AGE_BIN"] == age) &\
                                  (mut_total_count["TISSUE"] == tissue) &\
                                       (mut_total_count["BIN"] == bin_freq)]["TOTAL_MUT_COUNT"].unique()[0]

if check_flag == 0:
    #create the empty df with just the corresponding column fields
    output_empty_df =  pd.DataFrame({"STRAIN" : pd.Series(dtype='str'),
              "TISSUE" : pd.Series(dtype='str'),
              "AGE_BIN" : pd.Series(dtype='str'),
              "BIN" : pd.Series(dtype='str'),
              "SIM_RUN" : pd.Series(dtype='int'),
              "GENE" : pd.Series(dtype='str'),
              "ANNOTATION" : pd.Series(dtype='str'),
              "COUNT" : pd.Series(dtype='int'),
})
    #output the empty df 
    output_empty_filename = "output/indiv_sims/{}_{}_{}_{}_annotated_simulations_counts_per_gene.txt".format(strain, tissue, age, bin_freq)
    
    output_empty_df.to_csv(output_empty_filename, header = True, index = False, sep = "\t")
   
    #stop executing the program since we don't need to proceed
    sys.exit("Total count for this exp condition is 0, we will exit the program and save an empty file")

#load our annotation file of all possible variants
annotation_file = args.annotated_variants_file_path
annotations = pd.read_csv(annotation_file, sep = "\t")

#load our chrM sequence and create our dictionary of positions
chrM_file = args.fasta_file_path
#reading a fasta file in this makes our ref_seq a SeqRecord in biopython -- we treat the chrM ref seq as 1 record
ref_seq_record = SeqIO.read(chrM_file, "fasta")
ref_seq = ref_seq_record.seq.upper()

#merging our files to have the information we need in one position

all_conditions_info = pd.merge(mut_type_prop, mut_total_count, how = "left", on = ["STRAIN",\
                                                                                 "TISSUE",\
                                                                                 "AGE_BIN",\
                                                                                 "BIN"])

#load file needed to make the haplotype dictionary
haplotype_info_file = args.haplotype_file_path
haplotype_info = pd.read_csv(haplotype_info_file, sep = "\t")

#only keep SNVs
haplotype_snvs_only = haplotype_info[haplotype_info["VARIANT_TYPE"] == "SNV"]

#filter out haplotypes that are in the D-Loop and drop duplicates from our multiple experimental conditions 
haplotype_wo_dloop = haplotype_snvs_only[haplotype_snvs_only["START"] < 15422].drop_duplicates(subset = ["STRAIN", "START"], keep = "first")[["STRAIN","START", "REF", "ALT"]]

strains = list(haplotype_wo_dloop["STRAIN"].unique())

haplotype_dict = {}
for conplastic in strains:
    #take the position, ref, and alt alleles for all haploptypes for a given strain & place into one list of lists
    haplotypes_for_strain = haplotype_wo_dloop[haplotype_wo_dloop["STRAIN"] == conplastic]
    haplotypes_for_strain = haplotypes_for_strain[["START", "REF", "ALT"]].values.tolist()
    #add the list of lists as a value to the strain keys
    haplotype_dict[conplastic] = haplotypes_for_strain

#functions we will be using

#strain is defined via user input and haplotype_dict is created above
def get_ref_genome(strain, haplotype_dict):
   #we need to convert from a seq to a mutable object in order to mutate our sequence
    #recall ref_seq is a SeqRecord that was created using SeqIO.read on a fasta file of the chrM genome
    #our strain ref genome is the seq object associated with the SeqRecord
    mutable_strain_ref_genome = MutableSeq(ref_seq)

    haplotypes = haplotype_dict[strain]

    #we loop through each haplotype and mutate the ref sequence from mm10 accordingly
    for index in range(0, len(haplotypes)):
        site = haplotypes[index][0]
        allele = haplotypes[index][2]
        
        #a check to make sure that the site we are mutating is correct 
        if mutable_strain_ref_genome[site] != haplotypes[index][1]:
            print("There's an error in identifying the ref genome site") 
            
        #changing the allele at position in the ref genome     
        mutable_strain_ref_genome[site] = allele

    #convert back to a seq object in order to access all the biological methods 
    strain_ref_genome = Seq(mutable_strain_ref_genome)
	
    return strain_ref_genome

#this function takes as input a seq object, which is the object we mutated in the get_ref_genome
def get_position_dict(strain_ref_genome):
    #create a dictionary to hold the positions that our reference alleles are in
    #this will be used to choose which places in the genome are randomly mutated
    position_dict = {}

    nucleotides = ['A','C','G','T']

    #create the position_dict which the following format {"nucleotide": [lst of positions this nucleotide is located]}

    for base in nucleotides:
        position_dict[base] = [index for index in range(len(strain_ref_genome)) if strain_ref_genome[index] == base]

    return position_dict

#obtain the conplastic strain reference genome
if strain == "B6":
    strain_ref_genome = ref_seq
else:
    strain_ref_genome = get_ref_genome(strain, haplotype_dict)

#create the position_dict based on the conplastic strain ref genome
position_dict = get_position_dict(strain_ref_genome)

#subsetting our df with all the information for the coniditon we are processing at the moment
condition_info_df = all_conditions_info[(all_conditions_info["STRAIN"] == strain) &\
                                  (all_conditions_info["AGE_BIN"] == age) &\
                                  (all_conditions_info["TISSUE"] == tissue) &\
                                       (all_conditions_info["BIN"] == bin_freq)]

mut_count_array = condition_info_df["TOTAL_MUT_COUNT"].unique()

if len(mut_count_array) > 1:
    print("ERROR OCCURRED IN CALCULATING MUT_COUNT")

#getting integer out of the array 
mut_count = mut_count_array[0]

#creating a dictionary of our mutation type and heteroplasmy prop to call in in our multinomal sampling
mut_type_prop_dict = dict(zip(condition_info_df.MUT_TYPE, condition_info_df.PROP))

#make sure that round off error does not break multinomial draw
sum_mut_type_prop = sum(mut_type_prop_dict.values())

#normalizes our mutation type proportions to ensure that our proportions sum to 1
mut_type_prop_dict.update((key, (value/sum_mut_type_prop)) for key,value in mut_type_prop_dict.items())

#it's redundant I KNOOOOW but I want to make sure we keep track of the prop wrt mutation type
sim_mut_draws = np.random.multinomial(mut_count, [mut_type_prop_dict["A>C"], mut_type_prop_dict["A>G"], mut_type_prop_dict["A>T"],\
                                                      mut_type_prop_dict["C>A"], mut_type_prop_dict["C>G"], mut_type_prop_dict["C>T"],\
                                                      mut_type_prop_dict["G>A"], mut_type_prop_dict["G>C"], mut_type_prop_dict["G>T"],\
                                                      mut_type_prop_dict["T>A"], mut_type_prop_dict["T>C"], mut_type_prop_dict["T>G"]],\
                                      sim_num)
list_of_mutations = list(mut_type_prop_dict.keys())

#generate the mutations given in each simulation
simulated_variants_list = []

counter = 0

for sim_index in range(len(sim_mut_draws)):

    #for every possible reference allele create a dataframe that includes:
    #the mutation type, strain, age, tissue, bin_freq, and simulation
    A_ref = pd.DataFrame()
    A_ref["MUT_TYPE"] = pd.Series(list_of_mutations[0:3]).repeat(sim_mut_draws[sim_index][0:3])
    A_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][0:3])).values
    A_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][0:3])).values
    A_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][0:3])).values
    A_ref["BIN"] = pd.Series(bin_freq).repeat(sum(sim_mut_draws[sim_index][0:3])).values
    A_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][0:3])).values

    C_ref = pd.DataFrame()
    C_ref["MUT_TYPE"] = pd.Series(list_of_mutations[3:6]).repeat(sim_mut_draws[sim_index][3:6])
    C_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][3:6])).values
    C_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][3:6])).values
    C_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][3:6])).values
    C_ref["BIN"] = pd.Series(bin_freq).repeat(sum(sim_mut_draws[sim_index][3:6])).values
    C_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][3:6])).values

    G_ref = pd.DataFrame()
    G_ref["MUT_TYPE"] = pd.Series(list_of_mutations[6:9]).repeat(sim_mut_draws[sim_index][6:9])
    G_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][6:9])).values
    G_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][6:9])).values
    G_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][6:9])).values
    G_ref["BIN"] = pd.Series(bin_freq).repeat(sum(sim_mut_draws[sim_index][6:9])).values
    G_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][6:9])).values

    T_ref = pd.DataFrame()
    T_ref["MUT_TYPE"] = pd.Series(list_of_mutations[9:12]).repeat(sim_mut_draws[sim_index][9:12])
    T_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][9:12])).values
    T_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][9:12])).values
    T_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][9:12])).values
    T_ref["BIN"] = pd.Series(bin_freq).repeat(sum(sim_mut_draws[sim_index][9:12])).values
    T_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][9:12])).values

    #drawing uniformly with replacement from the position vectors

    A_ref_positions_mutated = [np.random.choice(position_dict["A"], \
                                                             sim_mut_draws[sim_index][draw]) \
                                                             for draw in range(0,3)]
    A_ref["POS"] = np.concatenate(A_ref_positions_mutated)

    C_ref_positions_mutated = [np.random.choice(position_dict["C"], \
                                                             sim_mut_draws[sim_index][draw]) \
                                                             for draw in range(3,6)]
    C_ref["POS"] = np.concatenate(C_ref_positions_mutated)

    G_ref_positions_mutated = [np.random.choice(position_dict["G"], \
                                                             sim_mut_draws[sim_index][draw]) \
                                                             for draw in range(6,9)]
    G_ref["POS"] = np.concatenate(G_ref_positions_mutated)

    T_ref_positions_mutated = [np.random.choice(position_dict["T"], \
                                                             sim_mut_draws[sim_index][draw]) \
                                                             for draw in range(9,12)]
    T_ref["POS"] = np.concatenate(T_ref_positions_mutated)

    #appending the above dataframes together -- this is much faster than the command commented out below for future reference
    simulated_variants_for_condition = pd.concat([A_ref, C_ref, G_ref, T_ref])

    simulated_variants_for_condition[["REF","ALT"]] = simulated_variants_for_condition["MUT_TYPE"].str.split(">", expand = True)

    simulated_variants_list.append(simulated_variants_for_condition)

    if counter % 1000 == 0:
        print(strain, age, tissue, bin_freq, sim_index)

    counter += 1

simulated_variants_df = pd.concat(simulated_variants_list)

filtering_noncoding_regions = simulated_variants_df[simulated_variants_df["POS"] >= 2750] 

annotated_simulated_variants = pd.merge(filtering_noncoding_regions, annotations, \
                                        how = "left", \
                                        on = ["STRAIN", "REF", "ALT", "POS"])[["STRAIN", "TISSUE", "AGE_BIN", "BIN",\
                                                                              "SIM_RUN", "POS", "REF", "ALT",\
                                                                              "MUT_TYPE", "ANNOTATION", "GENE",\
                                                                              "CODON_INDEX", "CODON_POS", "REF_AA", "ALT_AA"]]

protein_coding_annotated_simulated_variants = annotated_simulated_variants[annotated_simulated_variants["ANNOTATION"] != "non_coding_transcript_exon_variant"]

annotation_counts = protein_coding_annotated_simulated_variants[["STRAIN", "TISSUE", "AGE_BIN", "BIN",\
                                                  "SIM_RUN", "GENE", "MUT_TYPE", \
                                                  "ANNOTATION"]].groupby(["STRAIN", "TISSUE", "AGE_BIN", "BIN", "SIM_RUN", "GENE", "ANNOTATION"]).size()

annotation_counts = annotation_counts.reset_index()
annotation_counts.rename(columns = {0: "COUNT"}, inplace = True)

#creating the file name to output our simulations 
output_annotated_simulated_variants_per_gene = "output/indiv_sims/{}_{}_{}_{}_annotated_simulations_counts_per_gene.txt".format(strain, tissue, age, bin_freq)
output_raw_simulated_vars_w_annotations = "output/indiv_sims/{}_{}_{}_{}_raw_simulated_vars_w_annotations.txt".format(strain, tissue, age, bin_freq)

print("Exporting our files ...")
annotation_counts.to_csv(output_annotated_simulated_variants_per_gene, header = True, index = False, sep = "\t")
protein_coding_annotated_simulated_variants.to_csv(output_raw_simulated_vars_w_annotations, header = True, index = False, sep = "\t")
