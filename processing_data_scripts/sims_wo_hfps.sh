#!/bin/bash

python hNhS_simulator.py -strain "B6" -mut_count_file_path "files/mut_count_wo_hfps.txt" -mut_prop_file_path "files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -strain "AKR" -mut_count_file_path "files/mut_count_wo_hfps.txt" -mut_prop_file_path "files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -strain "ALR" -mut_count_file_path "files/mut_count_wo_hfps.txt" -mut_prop_file_path "files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -strain "F" -mut_count_file_path "files/mut_count_wo_hfps.txt" -mut_prop_file_path "files/mut_prop_wo_hfps.txt" -wo_hfp

python hNhS_simulator.py -strain "NZB" -mut_count_file_path "files/mut_count_wo_hfps.txt" -mut_prop_file_path "files/mut_prop_wo_hfps.txt" -wo_hfp
