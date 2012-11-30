Functional Tests to make sure Spats is trimming and mapping reads appropriately
-------------------------------------------------------------------------------
Requirements:

Bowtie 0.12.8
FASTX

File Description:
Reads_to_Map.txt - A list of constructed reads that should be mapped '(Keep)' if everything is working properly.
R1.fastq, R2.fastq - Reads in fastq format used by the test.
Targets.fa - Input targets file for the test.
Positive/ - Spats outputs for correctly trimmed and mapped. Everything listed as 'Keep' in Reads_to_Map should be mapped. See RRRY.sam and YYYR.sam for expected output if things are working correctly.
Negative/ - Output expected if adapter trimming and mapping is not working on the system. See RRRY.sam and YYYR.sam for expected output if things are NOT working correctly.
ApE/ - Ape files for constructing all the different read scenarios.

Example runtime outputs for both cases and the commands used.

Positive Test
=============

spats -p 2 --adapter-t AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter-b AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --trim-match 9 -o Positive Targets.fa RRRY YYYR R1.fastq R2.fastq

[Sun Nov 18 16:00:08 2012] Beginning Spats run (v0.0.1)
-----------------------------------------------
[Sun Nov 18 16:00:08 2012] Preparing output location Positive/
[Sun Nov 18 16:00:08 2012] Checking for Bowtie
	Bowtie version:		 0.12.8.0
[Sun Nov 18 16:00:08 2012] Relabeling reads in R1.fastq and R2.fastq
Kept 14 of 14 reads
[Sun Nov 18 16:00:08 2012] Indexing target sequences
[Sun Nov 18 16:00:09 2012] Trimming adapters from reads in Positive//NOMASK_1.fq
[Sun Nov 18 16:00:09 2012] Trimming adapters from reads in Positive//NOMASK_2.fq
[Sun Nov 18 16:00:09 2012] Rematching read pairs
[Sun Nov 18 16:00:09 2012] Relabeling reads in Positive/tmp/NOMASK_1.kept.fq and Positive/tmp/NOMASK_2.kept.fq
Kept 14 of 14 reads
[Sun Nov 18 16:00:09 2012] Mapping reads against targets with Bowtie
[Sun Nov 18 16:00:09 2012] Mapping reads against targets with Bowtie
[Sun Nov 18 16:00:09 2012] Building reactivity profiles
Processed 10 properly paired fragments, kept 5/5 (1.000000) treated, 5/5 (1.000000) untreated
-----------------------------------------------
Run complete [00:00:01 elapsed]

Negative Test
=============
spats -p 2 -o Negative Targets.fa RRRY YYYR R1.fastq R2.fastq
[Sun Nov 18 16:01:21 2012] Beginning Spats run (v0.0.1)
-----------------------------------------------
[Sun Nov 18 16:01:21 2012] Preparing output location Negative/
[Sun Nov 18 16:01:21 2012] Checking for Bowtie
	Bowtie version:		 0.12.8.0
[Sun Nov 18 16:01:21 2012] Relabeling reads in R1.fastq and R2.fastq
Kept 14 of 14 reads
[Sun Nov 18 16:01:21 2012] Indexing target sequences
[Sun Nov 18 16:01:21 2012] Relabeling reads in Negative//NOMASK_1.fq and Negative//NOMASK_2.fq
Kept 14 of 14 reads
[Sun Nov 18 16:01:21 2012] Mapping reads against targets with Bowtie
[Sun Nov 18 16:01:21 2012] Mapping reads against targets with Bowtie
[Sun Nov 18 16:01:21 2012] Building reactivity profiles
Processed 6 properly paired fragments, kept 3/3 (1.000000) treated, 3/3 (1.000000) untreated
-----------------------------------------------
Run complete [00:00:00 elapsed]



