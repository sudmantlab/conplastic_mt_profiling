library(tidyverse)
library(pryr)

outdir_files = "../files/"

print("Importing needed files ...")

possible_sites_per_gene_file = "../files/total_annotated_num_sites_per_gene.txt"
possible_site_per_gene = read.table(possible_sites_per_gene_file, header = TRUE, stringsAsFactors = FALSE)

simulations_file = "../files/all_sims_wo_dloop_sims.txt"
raw_simulations = read.table(simulations_file, header = TRUE, stringsAsFactors = FALSE)

print("Reformatting our raw simulations file ... ")

simulations_wide = raw_simulations %>%
	mutate(STRAIN = ifelse(STRAIN == "F", "FVB", STRAIN)) %>%
	pivot_wider(names_from = ANNOTATION, values_from = COUNT) %>%
	rename("SIM_NONSYN_MUT_COUNT" = "missense_variant", "SIM_SYN_MUT_COUNT" = "synonymous_variant")

#in the case that we did not simulate any syn or nonsyn mutations
simulations_wide[is.na(simulations_wide)] = 0

rm(raw_simulations)

print("Merging our sims data and our total number of possible position types in genes ...") 

merged_df = simulations_wide %>%
	left_join(possible_site_per_gene, by = c("STRAIN", "TISSUE", "AGE_BIN", "GENE")) %>%
	select(SIM_RUN, STRAIN, TISSUE, AGE_BIN, GENE, SIM_NONSYN_MUT_COUNT, TOTAL_NONSYN_SITE_COUNT, SIM_SYN_MUT_COUNT, TOTAL_SYN_SITE_COUNT)

print("Calculating our ratios ...")
sim_hnhs = merged_df %>%
	mutate(SIM_HN = SIM_NONSYN_MUT_COUNT/TOTAL_NONSYN_SITE_COUNT, SIM_HS = SIM_SYN_MUT_COUNT/TOTAL_SYN_SITE_COUNT) %>%
	mutate(SIM_HNHS = SIM_HN/SIM_HS)

print("Writing out our file ...")
mem_used()
write.table(sim_hnhs, file = paste(outdir_files,"sims_hNhS_ratios_per_gene_wo_dloop.txt", sep = ""), sep = "\t", quote = F, row.names = F)

print("Yay -- processing complete!")
