# spats
About Spats
===========

Spats is a package originally developed to support: 

- Julius B. Lucks, Stefanie A. Mortimer, Cole Trapnell, Shujun Luo, 
	Sharon Aviron, Gary P. Schroth, Lior Pachter, Jennifer A. Doudna, and Adam P. Arkin, 
	SHAPE-Seq: Multiplexed RNA Secondary and Tertiary Structure Determination 
        PNAS, 2011, doi: 10.1073/pnas.1106501108
- Sharon Aviran, Cole Trapnell, Julius B. Lucks, Stefanie A. Mortimer, Shujun Luo,
	Gary P. Schroth, Jennifer A. Doudna, Adam P. Arkin, and Lior Pachter, 
	Modeling and automation of SHAPE-Seq analysis 
        PNAS, 2011, doi: 10.1073/pnas.1106541108

and now supports:

- David Loughrey, Kyle E. Watters, Alexander H. Settle, and Julius B. Lucks
  SHAPE-Seq 2.0: Systematic Optimization and Extension of High-Throughput
  Chemical Probing of RNA Structure with Next Generation Sequencing
      Nucleic Acids Research, 2014 .
- Kyle E. Watters, Timothy R. Abbott, and Julius B. Lucks
  Simultaneous characterization of cellular RNA structure and function
  with in-cell SHAPE-Seq
      Nucleic Acids Research, 2015
- Kyle E. Watters, Angela M Yu, Eric J. Strobel, Alexander H. Settle, and Julius B. Lucks
  Characterizing RNA structures in vitro and in vivo with selective 2'-hydroxyl acylation
  analyzed by primer extension sequencing (SHAPE-Seq)
      Methods, 2016

- Watters KE, Strobel EJ, Yu AM, Lis JT, Lucks JB.
  Cotranscriptional folding of a riboswitch at nucleotide resolution
      Nature Structure and Molecular Biology 10.1038/nsmb.3316

and other papers from the SHAPE-Seq community.

The Spats package implements a read mapping and reactivity analysis pipeline for calculating
SHAPE-Seq reactivities from an input set of next-generation reads. The reactivity calculation
implemented is described in:

- Sharon Aviran, Julius B. Lucks, and Lior Pachter
  RNA Structure characterization from chemical mapping experiments.
  49th Allerton Conference, UIUC 2011. doi: 10.1109/Allerton.2011.6120379 .


For more information about the SHAPE-Seq experiment, and the application of spats for data analysis,
please see:

- Stephanie A. Mortimer, Cole Trapnell, Sharon Aviran, Lior Pachter, and Julius B. Lucks
  SHAPE-Seq: High Throughput RNA Structure Analysis
      Current Protocols in Chemical Biology, 2012, doi: 10.1002/9780470559277.ch120019

Spats Installation
==================

Check for latest instructions at: http://luckslab.github.io/spats/installation.html

* Requirements
  * Bowtie 0.12.8, must be available in your `PATH`. Get it at: http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/
  * cutadapt: python package (`pip install cutadapt`)

* Build from Source
  * `make -s`
  * add `bin` directory to your `PATH`


Spats Usage
===========
For detailed usage and examples, please see http://luckslab.github.io/spats/manual.html

The spats pipeline is run in two stages: adapter trimming and spats analysis. Adapter trimming is implemented in the adapter_trimmer code, while spats analysis is implemented in the spats toolkit. A typical workflow is as follows:

adapter_trimmer --A-b-sequence=<second_adapter_sequence> --A-t-sequence=<first_adapter_sequence> --read-len=35 R1.fastq R2.fastq RNA_targets.fa

spats --num-mismatches 0 -o Output RNA_targets.fa RRRY YYYR combined_R1.fastq combined_R2.fastq

where <second_adapter_sequence> is the sequence of the second adapter (present at the 3' end of cDNAs), <first_adapter_sequence> is the sequence of the first adapter (present at the 5' end of cDNAs), R1.fastq and R2.fastq are the fastq sequencing files from the sequencing run, and RNA_targets.fa is the FASTA formatted file containing the RNA sequences under study. Spats outputs text files (in Output) containing (+) and (-) fragment counts, beta's and theta’s, for each position in each RNA in the targets file.

Usage:
     adapter_trimmer [options]* <R1_seq.fastq> <R2_seq.fastq> <targets.fa>

     spats [options]* <targets.fa> <treated_handle> <untreated_handle> <R1_seq.fastq> <R2_seq.fastq>

     reactivities_split [options] <reactivities_file.out>

<targets.fa> is a fasta file containing the RNA sequences under study.
<R1_seq.fastq> is a fastq file containing "_1" reads from a paired-end Illumina run.
<R2_seq.fastq> is a fastq file containing "_2" reads from a paired-end Illumina run.
<treated_handle> and <untreated_handle> are IUPAC strings indicating the barcoding
scheme for the 1M7-treated and untreated pools.
<reactivities_file.out> is the output reactivities file from spats.

Given a paired end run (2x35bp) with target sequences (in DNA alphabet) in the file 
target_bcs.fa, "_1" reads in HCT20611_50bp_s_7_1_sequence.fq, and 
"_2" reads in HCT20611_50bp_s_7_2_sequence.fq, you can run the mapping pipeline like this:

adapter_trimmer --A-b-sequence=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --A-t-sequence=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC --read-len=35 HCT20611_50bp_s_7_1_sequence.fq HCT20611_50bp_s_7_2_sequence.fq target_bcs.fa

spats --num-mismatches 0 -o HCT20611_50bp.final RRRY YYYR combined_R1.fastq combined_R2.fastq

cd HCT20611_50bp.final
reactivities_split reactivities.out
