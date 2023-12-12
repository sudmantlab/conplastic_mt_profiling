#This script takes as input the vcf file for our simulations 
#To run this script run bash script running_per_gene_counts.sh

#Written by: Isabel M. Serrano on 08-16-2022

library(tidyverse)
library(argparse)

#Parsing arguments
parser = ArgumentParser()

parser$add_argument('--sims_file', help = "path to the simulated vcf files for the strain you want to process")
parser$add_argument('--strain', help = "strain we are processing info for")

xargs = parser$parse_args()

strain = xargs$strain

print(strain)
print("Importing vcf ...")

sims = read.table(xargs$sims_file, header=TRUE, stringsAsFactors = FALSE)

#Correcting FVB label if the sims file is the one corresponding to FVB
if (strain == "FVB"){
  sims = sims %>%
    mutate(STRAIN = "FVB")
}

print("Calculating count per gene for each sim run ...") 

count_df = sims %>% 
  select(SIM_RUN, STRAIN, TISSUE, AGE_BIN, GENE, ANNOTATION) %>%
  group_by(SIM_RUN, STRAIN, TISSUE, AGE_BIN, GENE, ANNOTATION) %>%
  summarise(COUNT = n())

rm(sims) 

print("Writing out file ...")
outdir_files = "output"

write.table(count_df, file = paste(outdir_files,"/",strain,"_count_per_gene.txt", sep = ""), sep = "\t", quote = F, row.names = F) 


print("Finished processing -- yay!")
