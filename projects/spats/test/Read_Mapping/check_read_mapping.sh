#!/usr/bin/bash
# Test cases for read mapping in spats
# Usage: bash check_read_mapping.sh
# Version: 0.1
# Created by Angela M Yu, March 29, 2016
# Copyright (c) 2016 Lucks Laboratory - all rights reserved.

export LD_LIBRARY_PATH=/fs/europa/g_jbl/Software/SHAPE-Seq/bin/lib:$LD_LIBRARY_PATH
export BOOST_ROOT=/fs/europa/g_jbl/Software/SHAPE-Seq/build/boost_1_49_0
export PYTHONPATH=$PYTHONPATH:/fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin
export PATH=/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/fs/home/amy35/bin:/opt/voyager/nbs/bin/:/usr/local/R/icse/bin:/fs/home/amy35/Projects/spats/build/spats-0.8.0/SHAPE-Seq/bin:/fs/home/amy35/Projects/spats/build/spats-0.8.0/SHAPE-Seq/bin/bin:/fs/europa/g_jbl/Software/SHAPE-Seq/bin:/fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin

output_name="SRP_All_Stops"
spats_output_dir="Output_test"
test_results_file="./test_results.txt"
sequence="ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCC"
adapter="CTGACTCGGGCACCAAGGAC"
rev_comp_adapter="GTCCTTGGTGCCCGAGTCAG"

# Remove previous test results and create target FASTA file
rm $test_results_file
printf ">${output_name}_${#sequence}nt\n$sequence$rev_comp_adapter\n" > ${output_name}.fa

# Iterate through all 16 binary permutations of length 4. This represents the pattern of reads mapped we expect for every 2 reads in their + and - channel.
for permutation in $(seq 0 15)
do
 python read_constructor.py --output $output_name --length 35 --sequence $sequence --adapter $adapter --permutation $permutation

 /fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin/adapter_trimmer_cutadapt.py --trim-match 8 --A-b-sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --A-t-sequence=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --read-len=35 ${output_name}_R1.fq ${output_name}_R2.fq ${output_name}.fa

 OUTPUT="$(ls -t | head -1)"

 /fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin/spats_cutadapt --all-RT-starts --num-mismatches 0 -o $spats_output_dir ${output_name}.fa RRRY YYYR ./$OUTPUT/combined_R1.fastq ./$OUTPUT/combined_R2.fastq

 python check_reactivities.py --input $spats_output_dir/reactivities.out --output $test_results_file --permutation $permutation --adapter $adapter --sequence $sequence

 rm -r $OUTPUT/ $spats_output_dir
 
done

# Output results and remove temporary files
echo -e "----------------------\nResults:"
cat $test_results_file
rm ${output_name}.fa ${output_name}_R1.fq ${output_name}_R2.fq
if ls ./targets_*/ 1> /dev/null 2>&1; then
 rm -r ./targets_*/
fi
