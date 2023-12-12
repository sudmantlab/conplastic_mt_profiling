#!/bin/bash
#this script concatenates all of our simulation files into one -- recall we had exporting issues due to memory limitations

head -1 output/indiv_sims/AKR_Brain_OLD_upper-5e-5_annotated_simulations_counts_per_gene.txt > output/final/combined_simulations_counts_per_gene.txt

for file in output/indiv_sims/*_annotated_simulations_counts_per_gene.txt
do
echo $file
tail -n +2 $file >> output/final/combined_simulations_counts_per_gene.txt 
done;

