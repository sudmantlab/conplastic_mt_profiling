#this script was modified on 07/26/2023 by IMS to run the simulator without the D-Loop. Modifications include: fasta file path change to a new fasta file without d-loop info, deleting filtering of the D-Loop region since we do that already, and a change to the output file name to include wo_dloop in file name

#modified on 07/31/2023 by IMS: the haplotype dict was edited such that we exclude the D-Loop variants
#include an if statement in our loop to skip the B6, Young, Liver iteration since we do not have data for this condition.
#added age as an input argument rather than a for loop due to memory issues in exporting sims
#added a sim_num argument to all for an easy variability in the number of simulations we want to run 
#added the fasta_file_path where the user can specify the fasta file they want to use with ease

import pandas as pd
import numpy as np
import io
import os
import argparse
from Bio import SeqIO

#adding user defined inputs for our simulator
parser = argparse.ArgumentParser()
parser.add_argument("-strain", help = "Conplastic Strain Name, IN ALL CAPS")
parser.add_argument("-age", help = "Age bin for samples, IN ALL CAPS")
parser.add_argument("-sim_num", help = "Number of sims you wish to generate")
parser.add_argument("-fasta_file_path", help = "File path to the corresponding fasta file")
parser.add_argument("-mut_count_file_path", help = "File path to the mut count file you are referring to")
parser.add_argument("-mut_prop_file_path", help = "File path to the mut prop file you are referring to")
parser.add_argument("-wo_hfp", action = "store_true", help = "Includes wo_hfp in output")
args = parser.parse_args()

#parsing our arguments to use in our simulator
strain = args.strain
age = args.age
sim_num = int(args.sim_num)

if args.wo_hfp:
    print("wo_hfps")
else:
    print("w_hfps")

#read in our mut_prop and mut_total_count files
mut_prop_file = args.mut_prop_file_path 
mut_total_count_file = args.mut_count_file_path 

mut_type_prop = pd.read_csv(mut_prop_file, sep = "\t")
mut_total_count = pd.read_csv(mut_total_count_file, sep = "\t")

#loading in our annotation file of all possible variants
annotation_file = "files/annotated_all_possible_variants.txt"
annotations = pd.read_csv(annotation_file, sep = "\t")

#read in our chrM sequence and create our dictionary of positions
chrM_file = args.fasta_file_path
ref_seq = SeqIO.read(chrM_file, "fasta")

#merging our files to have the information we need in one position

per_condition_info = pd.merge(mut_type_prop, mut_total_count, how = "left", on = ["STRAIN",\
                                                                                 "TISSUE",\
                                                                                 "AGE_BIN"])

#file needed to make the haplotype dictionary
haplotype_info_file = "files/haplotype_mutations.vcf"
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
    
def get_position_dict(strain_ref_genome):
    #create a dictionary to hold the positions that our reference alleles are in
    #this will be used to choose which places in the genome are randomly mutated
    position_dict = {}

    nucleotides = ['A','C','G','T']

    #create the position_dict which the following format {"nucleotide": [lst of positions this nucleotide is located]}

    for base in nucleotides:
        position_dict[base] = [index for index in range(len(strain_ref_genome)) if strain_ref_genome[index] == base]

    return position_dict

simulated_variants_list = []

#obtain the conplastic strain reference genome
if strain == "B6":
    strain_ref_genome = ref_seq
else:
    strain_ref_genome = get_ref_genome(strain, haplotype_dict)

#create the position_dict based on the conplastic strain ref genome
position_dict = get_position_dict(strain_ref_genome)

#strain and age were both passed in as user-input parameters
for tissue in ["Brain", "Heart", "Liver"]:

    if strain == "B6" and age == "YOUNG" and tissue == "Liver":
        continue

    condition_df = per_condition_info[(per_condition_info["STRAIN"] == strain) &\
                                      (per_condition_info["AGE_BIN"] == age) &\
                                      (per_condition_info["TISSUE"] == tissue)]

    mut_count = condition_df[(condition_df["AGE_BIN"] == age) \
                                   & (condition_df["TISSUE"] == tissue)]["HET_COUNT"].unique()

    if len(mut_count) > 1:
        print("ERROR OCCURRED IN CALCULATING MUT_COUNT")

    mut_count = mut_count[0]

    #creating a dictionary of our mutation type and heterplasmy prop to call in in our multinomal sampling
    mut_type_prop_dict = dict(zip(condition_df.MUTATION_TYPE, condition_df.HET_PROP))

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

    counter = 0

    for sim_index in range(len(sim_mut_draws)):

        #for every possible reference allele create a dataframe that includes:
        #the mutation type, strain, age, tissue, simulation
        A_ref = pd.DataFrame()
        A_ref["MUT_TYPE"] = pd.Series(list_of_mutations[0:3]).repeat(sim_mut_draws[sim_index][0:3])
        A_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][0:3])).values
        A_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][0:3])).values
        A_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][0:3])).values
        A_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][0:3])).values

        C_ref = pd.DataFrame()
        C_ref["MUT_TYPE"] = pd.Series(list_of_mutations[3:6]).repeat(sim_mut_draws[sim_index][3:6])
        C_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][3:6])).values
        C_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][3:6])).values
        C_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][3:6])).values
        C_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][3:6])).values

        G_ref = pd.DataFrame()
        G_ref["MUT_TYPE"] = pd.Series(list_of_mutations[6:9]).repeat(sim_mut_draws[sim_index][6:9])
        G_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][6:9])).values
        G_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][6:9])).values
        G_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][6:9])).values
        G_ref["SIM_RUN"] = pd.Series(sim_index).repeat(sum(sim_mut_draws[sim_index][6:9])).values

        T_ref = pd.DataFrame()
        T_ref["MUT_TYPE"] = pd.Series(list_of_mutations[9:12]).repeat(sim_mut_draws[sim_index][9:12])
        T_ref["STRAIN"] = pd.Series(strain).repeat(sum(sim_mut_draws[sim_index][9:12])).values
        T_ref["AGE_BIN"] = pd.Series(age).repeat(sum(sim_mut_draws[sim_index][9:12])).values
        T_ref["TISSUE"] = pd.Series(tissue).repeat(sum(sim_mut_draws[sim_index][9:12])).values
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

        #joined_variants = A_ref.append(C_ref).append(G_ref).append(T_ref).reset_index(drop = True)

        simulated_variants_for_condition[["REF","ALT"]] = simulated_variants_for_condition["MUT_TYPE"].str.split(">", expand = True)

        simulated_variants_list.append(simulated_variants_for_condition)

        if counter % 1000 == 0:
            print(strain, age, tissue, sim_index)

        counter += 1

simulated_variants_df = pd.concat(simulated_variants_list)

filtering_noncoding_regions = simulated_variants_df[simulated_variants_df["POS"] >= 2750] 

annotated_simulated_variants = pd.merge(filtering_noncoding_regions, annotations, \
                                        how = "left", \
                                        on = ["STRAIN", "REF", "ALT", "POS"])[["STRAIN", "TISSUE", "AGE_BIN",\
                                                                              "SIM_RUN", "POS", "REF", "ALT",\
                                                                              "MUT_TYPE", "ANNOTATION", "GENE",\
                                                                              "CODON_INDEX", "CODON_POS", "REF_AA", "ALT_AA"]]

protein_coding_annotated_simulated_variants = annotated_simulated_variants[annotated_simulated_variants["ANNOTATION"] != "non_coding_transcript_exon_variant"]

annotation_counts = protein_coding_annotated_simulated_variants[["STRAIN", "TISSUE", "AGE_BIN",\
                                                  "SIM_RUN", "GENE", "MUT_TYPE", \
                                                  "ANNOTATION"]].groupby(["STRAIN", "TISSUE", "AGE_BIN", "SIM_RUN", "GENE", "ANNOTATION"]).size()

annotation_counts = annotation_counts.reset_index()
annotation_counts.rename(columns = {0: "COUNT"}, inplace = True)

if args.wo_hfp:
    print("wo hfp")
    output_annotated_simulated_variants = "output/{}_{}_wo_dloop_wo_hfp_sim_annotated_counts.txt".format(strain, age)
else:
    print("has hfp")
    output_annotated_simulated_variants = "output/{}_{}_wo_dloop_sim_annotation_counts.txt".format(strain, age)

annotation_counts.to_csv(output_annotated_simulated_variants, header = True, index = False, sep = "\t")
