#Updated by IMS 08-22-2023

import pandas as pd 
import numpy as np
#libraries in base python
import io
import os
import math
import argparse
#libraries we need to import from BioPython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

#parsing the arguments (mainly files) we will be using in this analysis 
#adding user defined inputs for our simulator
parser = argparse.ArgumentParser()
#parameters defining our experimental unit
parser.add_argument("-vcf_file_path", help = "File path to the vcf we want to annotate")
parser.add_argument("-fasta_file_path", help = "File path to the corresponding fasta file")
parser.add_argument("-haplotype_file_path", help = "File path to a file containing information about mouse haplotypes that distinguish it from mm10")
parser.add_argument("-mt_coords_file_path", help = "File path to a coordinates file for the mm10 chrM")
parser.add_argument("-output_file", help = "Name of output file you'd like to designate")

args = parser.parse_args()

##Import, filter, and format all of the files we will be working from##

#upload file we want to annotate
vcf_file = args.vcf_file_path
vcf = pd.read_csv(vcf_file,  sep = "\t")

#upload fasta file
chrM_file = args.fasta_file_path
ref_seq_record = SeqIO.read(chrM_file, "fasta")
#ref_seq is now a seq object that we made sure was all uppercase for downstream string comparison
ref_seq = ref_seq_record.seq.upper()

#file needed to make the haplotype dictionary 
haplotype_info_file = args.haplotype_file_path
haplotype_info = pd.read_csv(haplotype_info_file,  sep = "\t")

#filter out indels from the haplotype info (won't be included in our dictionary to mutate)
#only keep SNVs
haplotype_snvs_only = haplotype_info[haplotype_info["VARIANT_TYPE"] == "SNV"]

#file that contains information for the genes in the mtGenome
mt_dna_coords_file = args.mt_coords_file_path
mt_dna_coords = pd.read_csv(mt_dna_coords_file,  sep = "\t", header = None)

#Formatting the mt_dna_coords file 
#1) 0-indexing our mtDNA coordinates -- since Python 0 indexes 
#2) Recall that list slicing has an exclusive end so we don't need to update the end coordinates
#3) Adding the D-Loop to the coordinates

mt_dna_coords.rename(columns = {0: "GENE", 1: "START", 2:"END"}, inplace = True)
mt_dna_coords["START"] = mt_dna_coords["START"] - 1

#add the d-loop to the coordinates
mt_dna_coords.loc[len(mt_dna_coords.index)] = ['D-Loop',15422 , 16299]

#annotate the missing bp to avoid out of index errors
missing_bp = pd.DataFrame({"GENE": ["ND1-Ti", "Tq-Tm", "Tw-Ta", "Ta-Tn","Tc-Ty", "Ty-CO1", "Ts1-Td", "Td-CO2","CO2-Tk","Tk-ATP8", "ND3-Tr", "Tr-ND4L", "ND4-Th", "Te-CytB"],\
                           "START": [3706, 3842, 5016, 5086, 5257, 5326, 6938, 7011, 7696, 7764, 9806, 9875, 11544, 14139],\
                           "END": [3707, 3844, 5017, 5088, 5259, 5327, 6941, 7012, 7699, 7765, 9807, 9876, 11545, 14144]})

mt_dna_coords = pd.concat([mt_dna_coords, missing_bp], ignore_index = True)

#create a column that flags protein coding regions
protein_coding = ["mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", \
         "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",\
        "mt-Nd5", "mt-Nd6", "mt-Cytb"]

mt_dna_coords["CODING"] = np.where(mt_dna_coords["GENE"].isin(protein_coding), "PROTEIN", "OTHER")

##create a dictionary with the following structure {"strain":[[site, ref, alt], [site, ref, alt]} -- builds the dictionary as we want##
haplotype_wo_dloop = haplotype_snvs_only[haplotype_snvs_only["START"] < 15422].drop_duplicates(subset =
 ["STRAIN", "START"], keep = "first")[["STRAIN","START", "REF", "ALT"]]

strains = list(haplotype_wo_dloop["STRAIN"].unique())

haplotype_dict = {}
for conplastic in strains:
    #take the position, ref, and alt alleles for all haploptypes for a given strain & place into one list of lists
    haplotypes_for_strain = haplotype_wo_dloop[haplotype_wo_dloop["STRAIN"] == conplastic]
    haplotypes_for_strain = haplotypes_for_strain[["START", "REF", "ALT"]].values.tolist()
    #add the list of lists as a value to the strain keys
    haplotype_dict[conplastic] = haplotypes_for_strain

##create the functions that we will be using to annotate our variants##
##get_ref_genome -- input: string of strain name; output: a seq object of the strain reference genome 
##get_gene_info -- input: mutation_pos (0-index) and the mt_dna_coords dataframe; output: a list of lists named region with the structure [[gene_name, gene_start, gene_end, gene_type]]
##get_codon_pos -- input: gene_start and mutation_pos; output: the position in the codon that the mutation is in 0, 1, or 2
##get_annotation -- input: strain_ref_genome, gene, gene_start, gene_end, mutation_pos, alt_allele; output: annotation, codon_index, codon_pos, ref_AA, alt_AA
##annotation_process -- input: sample, mutation_pos, ref_allele, alt_allele, strain_ref_genome, mt_dna_coords; output: appending to the results lists the following information sample, mutation_pos, ref_allele, alt_allele, annotation, gene, codon_index, codon_pos, ref_AA, alt_AA

def get_ref_genome(strain, haplotype_dict):
   #we need to convert from a seq to a mutable object in order to mutate our sequence
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

def get_gene_info(mutation_pos, mt_dna_coords):

    #identify the gene name, start, end, and coding type needed to annotate this position
    #in the case that there is more than one mtdna region that this mutation falls into, region will be a nested list
    #of length 2
    region = mt_dna_coords[(mt_dna_coords["START"] <= mutation_pos) \
                           & (mt_dna_coords["END"] > mutation_pos)].values.tolist()
    return region

def get_codon_position(mutation_pos, gene_start):
    codon_pos = (mutation_pos - gene_start)% 3
    return codon_pos

def get_annotation(strain_ref_genome, gene, gene_start, gene_end, mutation_pos, alt_allele):
    
    #obtain the reference gene sequence
    ref_gene = strain_ref_genome[gene_start:gene_end]
    
    #reverse complement ND6 -- the only gene on the negative strand all other annotations are 
    #with respect to the positive strand and the ref genome gives information with respect to the positive strand
    #note: reverse complimenting the ref_gene will give us the complement ref allele we need in the mutation position
    if gene == "mt-Nd6":
        ref_gene = ref_gene.reverse_complement()
        
    #convert to the amino acid sequence & convert to a string
    ref_amino_acid_seq = str(ref_gene.translate(table="Vertebrate Mitochondrial"))
    
    #mutate the reference genome to have the mutation we want to annotate
    #strain_ref_genome is a Seq object
    mutant_genome = MutableSeq(strain_ref_genome)
    mutant_genome[mutation_pos] = alt_allele
    #change the mutant genome back to a sequence object 
    mutant_genome_seq_obj = Seq(mutant_genome)
    #grab the mutatant gene
    mutant_gene = mutant_genome_seq_obj[gene_start:gene_end]
    
    #in the case that our mutatant gene is ND6 we need to reverse complement the gene
    #note: reverse complementing the mutant gene as a whole will complement the alternate allele 
    if gene == "mt-Nd6":
        mutant_gene = mutant_gene.reverse_complement()
    
    #translate the gene
    mutant_amino_acid_seq = str(mutant_gene.translate(table="Vertebrate Mitochondrial"))
    
    #now let's check if the strings are the same (synonymous variant) or if they differ and where (missense variant)
    mismatch = [i for i in range(0, len(mutant_amino_acid_seq)) \
                if mutant_amino_acid_seq[i] != ref_amino_acid_seq[i]]
    
    if mismatch:
        if len(mismatch) > 1:
            print("ERROR OCCURED - WE DON'T RETURN TO REF GENOME CORRECTLY")
            
        annotation = "missense_variant"
        codon_index = mismatch[0]
        ref_AA = ref_amino_acid_seq[codon_index]
        alt_AA = mutant_amino_acid_seq[codon_index]  
        codon_pos = get_codon_position(mutation_pos, gene_start)
        
        #we need to take into account that the end is now the beginning for ND6 and that we have end exclusion in slicing
        if gene == "mt-Nd6":
            codon_pos = ((gene_end - 1) - mutation_pos) % 3
    else:
        
        annotation = "synonymous_variant"
        codon_index = math.floor((mutation_pos - gene_start)/3)
        
        if (codon_index == len(ref_amino_acid_seq)) & (len(ref_gene) % 3):
            #this would be true for 3 transcripts COXIII, ND4, and CYTB
            codon_index = "LAST_CODON"
            codon_pos = "BEGIN_STOP_CODON"
            ref_AA = "*"
            alt_AA = "*"
            return annotation, codon_index, codon_pos, ref_AA, alt_AA
            
        codon_pos = get_codon_position(mutation_pos, gene_start)
        
        #we need to account for ND6 being "flipped around"; since we compare the amino acid sequences to get the 
        #codon index in missense mutations that case already accounts for reversing the strand
        if gene == "mt-Nd6":
            codon_index = math.floor(((gene_end - 1 ) - mutation_pos)/3)
            codon_pos = ((gene_end - 1) - mutation_pos) % 3
            
        ref_AA = ref_amino_acid_seq[codon_index]
        alt_AA = mutant_amino_acid_seq[codon_index]
        

    return annotation, codon_index, codon_pos, ref_AA, alt_AA

def annotating_process(strain, mutation_pos, ref_allele, alt_allele, strain_ref_genome, mt_dna_coords):
    
    #retrieve all the gene information for which the variant exists
    #check to see if the variant is in an overlapping region
    region = get_gene_info(mutation_pos, mt_dna_coords)
    
    if len(region) > 1:
        for gene_index in range(0, len(region)):
            gene, gene_start, gene_end, gene_type = region[gene_index][0],\
            region[gene_index][1], region[gene_index][2], region[gene_index][3]
            
            if gene_type == "OTHER": 
                annotation = "non_coding_transcript_exon_variant"
                codon_index = "NA"
                codon_pos = "NA"
                ref_AA = "NA"
                alt_AA = "NA"
            else: 
                annotation, codon_index, codon_pos, ref_AA, alt_AA = get_annotation(strain_ref_genome, gene,\
                                                                         gene_start, gene_end, mutation_pos, alt_allele)

            results.append([strain, mutation_pos, ref_allele, alt_allele, annotation, gene, codon_index,\
                           codon_pos, ref_AA, alt_AA])
    else: 
        gene, gene_start, gene_end, gene_type = region[0][0], region[0][1], region[0][2], region[0][3]
            
        if gene_type == "OTHER": 
            annotation = "non_coding_transcript_exon_variant"
            codon_index = "NA"
            codon_pos = "NA"
            ref_AA = "NA"
            alt_AA = "NA"
        else: 
            annotation, codon_index, codon_pos, ref_AA, alt_AA = get_annotation(strain_ref_genome, gene,\
                                                                     gene_start, gene_end, mutation_pos, alt_allele)

            #note we save the recording_pos which is 1-indexed and run the annotation on the 0-indexed location
            #of the variant
        results.append([strain, mutation_pos, ref_allele, alt_allele, annotation, gene, codon_index,\
                       codon_pos, ref_AA, alt_AA])

#Now, to start the annotation process we call all of our above functions
strains_to_check = sorted(vcf["STRAIN"].unique())

#This is the main body of the code that calls and applies our functions for each somatic snv
results = []

for strain in strains_to_check:
    
    #obtain the conplastic strain reference genome
    if strain == "B6":
        strain_ref_genome = ref_seq
    else:
        strain_ref_genome = get_ref_genome(strain, haplotype_dict)
    
    #subset the data to find all positions we need to check in the given strain
    subset_data = vcf[vcf["STRAIN"] == strain]
    
    subset_data.apply(lambda row: annotating_process(row["STRAIN"], row["POS"], row["REF"], row["ALT"], strain_ref_genome, mt_dna_coords), axis = 1)

#convert results into a dataframe for export 
annotations = pd.DataFrame(results, columns = ["STRAIN", "POS", "REF", "ALT", "ANNOTATION", "GENE", \
                                            "CODON_INDEX", "CODON_POS", "REF_AA", "ALT_AA"])

output_file = args.output_file
annotations.to_csv(output_file, header=True, index=False, sep = "\t")
