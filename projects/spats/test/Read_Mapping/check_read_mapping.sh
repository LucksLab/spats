#!/usr/bin/bash


export LD_LIBRARY_PATH=/fs/europa/g_jbl/Software/SHAPE-Seq/bin/lib:$LD_LIBRARY_PATH
export BOOST_ROOT=/fs/europa/g_jbl/Software/SHAPE-Seq/build/boost_1_49_0
export PYTHONPATH=$PYTHONPATH:/fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin
export PATH=/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/fs/home/amy35/bin:/opt/voyager/nbs/bin/:/usr/local/R/icse/bin:/fs/home/amy35/Projects/spats/build/spats-0.8.0/SHAPE-Seq/bin:/fs/home/amy35/Projects/spats/build/spats-0.8.0/SHAPE-Seq/bin/bin:/fs/europa/g_jbl/Software/SHAPE-Seq/bin:/fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin

output_name="SRP_All_Stops"
spats_output_dir="Output_test"
sequence="ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCC"
adapter="CTGACTCGGGCACCAAGGAC"
rev_comp_adapter="GTCCTTGGTGCCCGAGTCAG"

for permutation in $(seq 0 15)
do
 python read_constructor.py --output $output_name --length 35 --sequence $sequence --adapter $adapter --permutation $permutation
 printf ">${output_name}_${#sequence}nt\n$sequence$rev_comp_adapter\n" > ${output_name}.fa

 /fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin/adapter_trimmer_cutadapt.py --trim-match 8 --A-b-sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --A-t-sequence=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --read-len=35 ${output_name}_R1.fq ${output_name}_R2.fq ${output_name}.fa

 OUTPUT="$(ls -t | head -1)"

 /fs/europa/g_jbl/Software/SHAPE-Seq/bin/bin/spats_cutadapt --all-RT-starts --num-mismatches 0 -o $spats_output_dir ${output_name}.fa RRRY YYYR ./$OUTPUT/combined_R1.fastq ./$OUTPUT/combined_R2.fastq

 python check_reactivities.py --input $spats_output_dir/reactivities.out --permutation $permutation --adapter $adapter --sequence $sequence

 rm -r $OUTPUT/ $spats_output_dir
 
done

rm ${output_name}.fa ${output_name}_R1.fq ${output_name}_R2.fq
