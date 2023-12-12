#!/bin/bash
#this script concatenates all of our simulation files into one -- recall we had exporting issues due to memory limitations

head -1 output/AKR_YOUNG_wo_dloop_wo_hfp_sim_annotated_counts.txt > output/all_sims_wo_dloop_sims.txt

for file in output/*_wo_dloop_wo_hfp_sim_annotated_counts.txt
do
echo $file
tail -n +2 $file >> output/all_sims_wo_dloop_sims.txt
done;

