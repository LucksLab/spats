#!/usr/bin/bash
# Test case for read mapping in spats
# Usage: bash check_read_mapping.sh
# Version: 0.1
# Created by Angela M Yu, March 29, 2016
# Copyright (c) 2016 Lucks Laboratory - all rights reserved.

source /fs/europa/g_jbl/Software/Test/test_bash_profile

output_name="SRP_All_Stops"
spats_output_dir="Output_test"
test_results_file="./test_results.txt"
sequence="ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCC"
linker="CTGACTCGGGCACCAAGGAC"
rev_comp_linker="GTCCTTGGTGCCCGAGTCAG"

# Remove previous test results and create target FASTA file
rm $test_results_file
printf ">${output_name}_${#sequence}nt\n$sequence$rev_comp_linker\n" > ${output_name}.fa

python read_constructor.py --output $output_name --length 35 --sequence $sequence --linker $linker

adapter_trimmer --trim-match 8 --A-b-sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --A-t-sequence=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --read-len=35 ${output_name}_R1.fq ${output_name}_R2.fq ${output_name}.fa

OUTPUT="$(ls -t | head -1)"

spats --all-RT-starts --num-mismatches 0 -o $spats_output_dir ${output_name}.fa RRRY YYYR ./$OUTPUT/combined_R1.fastq ./$OUTPUT/combined_R2.fastq

python check_reactivities.py --input $spats_output_dir/reactivities.out --output $test_results_file --linker $linker --sequence $sequence

#Print results and remove combined_reads_*/
echo -e "\n-----------------------------------------------\nResults:"
cat test_results.txt
rm -r $OUTPUT/
