Read Mapping Test Suite
-----------------------

Tests read mapping by spats pipeline.

Files:

read_constructor.py - takes input sequences and constructs properly formatted (with Illumina sequences etc.) fastq reads for submission to spats.
test_read_mapping.sh - driver code to run the test
check_reactivities.py - code to check test. Used by test_read_mapping.sh

Running:

Note: Path must be set appropriately so that spats and adapter_trimmer can be reached.

bash test_read_mapping.sh

Expected Output:

[Tue Mar 29 21:09:22 2016] Beginning Adapter_trimmer reads processing (vTEST0.2.0)
-----------------------------------------------------
[Tue Mar 29 21:09:22 2016] Building Bowtie targets
[Tue Mar 29 21:09:23 2016] Determining optimal trim_match and minimum read length
[Tue Mar 29 21:09:23 2016] Automatically dropping 1 nt from end of RNAs (from 3' end analysis)
[Tue Mar 29 21:09:23 2016] Trim-match manually set to 8
[Tue Mar 29 21:09:23 2016] Filtering reads and storing Case I reads
[Tue Mar 29 21:09:23 2016] Clipping Case II reads
[Tue Mar 29 21:09:23 2016] Excecuting trim_search on remaining Case II reads
[Tue Mar 29 21:09:23 2016] Finding reads that are 30 nt long
[Tue Mar 29 21:09:24 2016] Finding reads that are 29 nt long
[Tue Mar 29 21:09:24 2016] Finding reads that are 28 nt long
[Tue Mar 29 21:09:24 2016] Finding reads that are 27 nt long
[Tue Mar 29 21:09:24 2016] Finding reads that are 26 nt long
[Tue Mar 29 21:09:24 2016] Finding reads that are 25 nt long
[Tue Mar 29 21:09:25 2016] Finding reads that are 24 nt long
[Tue Mar 29 21:09:25 2016] Combining all processed reads
-----------------------------------------------
Run complete [00:00:02 elapsed]

/fs/europa/g_jbl/Software/Test/Spats_Mapping_Test/combined_reads_LHJRS3PPZL/

[Tue Mar 29 21:09:25 2016] Beginning Spats run (vTEST1.0.0)
-----------------------------------------------
[Tue Mar 29 21:09:25 2016] Preparing output location Output_test/
[Tue Mar 29 21:09:25 2016] Checking for Bowtie
	Bowtie version:		 0.12.8.0
[Tue Mar 29 21:09:25 2016] Relabeling reads in ./combined_reads_LHJRS3PPZL/combined_R1.fastq and ./combined_reads_LHJRS3PPZL/combined_R2.fastq
Kept 14042 of 14042 reads
[Tue Mar 29 21:09:25 2016] Indexing target sequences
[Tue Mar 29 21:09:25 2016] Relabeling reads in Output_test//NOMASK_1.fq and Output_test//NOMASK_2.fq
Kept 14042 of 14042 reads
[Tue Mar 29 21:09:25 2016] Mapping reads against targets with Bowtie
[Tue Mar 29 21:09:26 2016] Mapping reads against targets with Bowtie
[Tue Mar 29 21:09:26 2016] Building reactivity profiles
Processed 14042 properly paired fragments, kept 7021/7021 (1.000000 %) treated, 7021/7021 (1.000000 %)  untreated
-----------------------------------------------
Run complete [00:00:00 elapsed]

-----------------------------------------------
Results:
OK - 272 read positions out of 272 expected, 272 correct, 0 incorrect