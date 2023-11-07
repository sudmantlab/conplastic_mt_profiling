#!/bin/bash

python hNhS_simulator.py -sim_num 10000 -strain "B6" -mut_count_file_path "../files/mut_count_wo_hfps.txt" -mut_prop_file_path "../files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -sim_num 10000 -strain "AKR" -mut_count_file_path "../files/mut_count_wo_hfps.txt" -mut_prop_file_path "../files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -sim_num 10000 -strain "ALR" -mut_count_file_path "../files/mut_count_wo_hfps.txt" -mut_prop_file_path "../files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -sim_num 10000 -strain "F" -mut_count_file_path "../files/mut_count_wo_hfps.txt" -mut_prop_file_path "../files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -sim_num 10000 -strain "NZB" -mut_count_file_path "../files/mut_count_wo_hfps.txt" -mut_prop_file_path "../files/mut_prop_wo_hfps.txt" -wo_hfp
