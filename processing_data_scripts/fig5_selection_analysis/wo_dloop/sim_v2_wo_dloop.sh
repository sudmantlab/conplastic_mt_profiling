#!/bin/bash

python hNhS_simulator_v2.py -strain "B6" -age "OLD" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "B6" -age "YOUNG" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "AKR" -age "OLD" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "AKR" -age "YOUNG" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "ALR" -age "OLD" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "ALR" -age "YOUNG" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "F" -age "OLD" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "F" -age "YOUNG" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "NZB" -age "OLD" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp

python hNhS_simulator_v2.py -strain "NZB" -age "YOUNG" -sim_num 10000 -fasta_file_path "files/chrM_wo_dloop.fa" -mut_count_file_path "files/total_mut_count_wo_dloop.txt" -mut_prop_file_path "files/mut_type_props_nondloop.txt" -wo_hfp
