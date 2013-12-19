# Test suite for adapter_trimmer.py
# File Descriptions
## test_suite.py - this file
## Reads Files:
### Reads files generated from the tRNA sequence (tRNA_full_length_with_adapters.ape)
### Reads are named according to the code: tRNA-XX-[Y/N]-[Y/N],
### where XX is how long the RNA fragment was, the first [Y/N] is whether the read should be 
### kept after the trim_match algorithm, and the second [Y/N] is whether the read should be 
### kept after the final trimming is performed by the adapter_trimmer algorithm
### Files:
#### R1_35bp_test.fq - fastq file containing R1 reads for tRNA test
#### R2_35bp_test.fq - fastq file containing R2 reads for tRNA test
#### Test_Correct/ - Fastq files representing a successful algorithm - for test comparisons
##### R1_35bp_test_trim_match_correct.fq - fastq file containing R1 reads for tRNA test that should be kept for trim_match positive test 
##### R2_35bp_test_trim_match_correct.fq - fastq file containing R2 reads for tRNA test that should be kept for trim_match positive test
##### R1_35bp_test_correct.fq - fastq file containing R1 reads for tRNA test that should be kept for full positive test 
##### R2_35bp_test_correct.fq - fastq file containing R2 reads for tRNA test that should be kept for full positive test
## targets_test.fa - targets file used for the test runs
## tRNA_full_length_with_adapters.ape - ape file containing tRNA sequence and adapters used to construct the test

# Running the test
## python test_suite.py
### Success
#### Trim_Match Test Passed
##### Full Test Passed
### Failure
#### Tests will report failed and test files will not be cleaned up. sdiff results will be reported for debugging

import os
import subprocess
import adapter_trimmer as at

R1_test_file = "R1_35bp_test.fq"
R2_test_file = "R2_35bp_test.fq"

class TestCases(object):
    """storing relevant info about these tests"""
    def __init__(self):
        super(TestCases, self).__init__()
        self.files = []
        self.directories = []
        
    def setup(self):
        """setting up the test"""
        # build bowtie index
        if not os.path.isdir("targets_index"):
            os.system("mkdir targets_index")
            self.directories.append("targets_index")
        os.system("bowtie-build -q targets_test.fa targets_index/targets")
    
        # build output directory
        if not os.path.isdir("test_output"):
            os.system("mkdir test_output")
            self.directories.append("test_output")
    
        # format reads to have perfect quality scores
        input_R1 = at.quality_mod(R1_test_file,35)
        input_R2 = at.quality_mod(R2_test_file,35)
        
        self.files.append(input_R1)
        self.files.append(input_R2)

    def cleanup(self,other_files=[],other_dirs=[]):
        """deleting temp files and directories"""
        for f in self.files:
            try:
                os.remove(f)
            except:
                pass
        for f in other_files:
            try:
                os.remove(f)
            except:
                pass
        for d in self.directories:
            try:
                os.system("rm -r %s" % d)
            except:
                pass
        for d in other_dirs:
            try:
                os.system("rm -r %s" % d)
            except:
                pass
            
def main():
    """Running the main tests"""
    
    # Trim Match Test
    trim_match_test = TestCases()
    trim_match_test.setup()
    # Test trim_search algorithm
    # Test case assuming 35 bp inserts and trim_match=5
    read_len = 35
    trim_match = 5
    trim_len = read_len-trim_match
    #for base in range(read_len,trim_len,-1): # wrong - picks up @tRNA-31-p-N-Y
    for base in range(read_len-1,trim_len,-1):
        at.trim_search(base,4,"test_output/","targets_index/",trim_match_test.files[0],trim_match_test.files[1])
    
    # Confirm Test - concatenate files in reverse order to put together a list of reads that kept
    os.system("cat `ls -r test_output/bowtie_results3*_1_rc.fq` > testR1")
    os.system("cat `ls -r test_output/bowtie_results3*_2.fq` > testR2")
    
    # Grab string of sdiff between generated test and expected results
    test_R1_check = subprocess.Popen(["sdiff","-s","testR1", "Test_Correct/R1_35bp_test_trim_match_correct.fq"],
                                     stdout=subprocess.PIPE).communicate()[0]
    test_R2_check = subprocess.Popen(["sdiff","-s","testR2", "Test_Correct/R2_35bp_test_trim_match_correct.fq"],
                                     stdout=subprocess.PIPE).communicate()[0]

    if test_R1_check == "" and test_R2_check == "": # Should see no differences between results and correct files
        print "Trim_Match Test Passed"
        trim_match_test.cleanup(["testR1","testR2"])
    else:
        print "Trim_Match Test Failed"
        print "sdiff -s testR1 Test_Correct/R1_35bp_test_trim_match_correct.fq"
        print test_R1_check
        print "sdiff -s testR2 Test_Correct/R2_35bp_test_trim_match_correct.fq"
        print test_R2_check
    
    
    # Full Test
    full_test = TestCases()
    full_test.setup()
    
    trim_match = 5
    read_len = 35
    A_b_sequence = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
    A_t_sequence = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    min_read_len = 1
    final_dir = "test_output/"
    trim_auto = False
    min_read_auto = False
    max_handle_len = 4
    input_R1 = "R1_35bp_test.fq"
    input_R2 = "R2_35bp_test.fq"
    input_targets = "targets_test.fa"
    at.full_trim(trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len,input_R1,input_R2,input_targets)
    
    # Grab string of sdiff between generated test and expected results
    test_R1_check = subprocess.Popen(["sdiff","-s",final_dir+"combined_R1.fastq", "Test_Correct/R1_35bp_test_correct.fq"],
                                     stdout=subprocess.PIPE).communicate()[0]
    test_R2_check = subprocess.Popen(["sdiff","-s",final_dir+"combined_R2.fastq", "Test_Correct/R2_35bp_test_correct.fq"],
                                     stdout=subprocess.PIPE).communicate()[0]
    if test_R1_check == "" and test_R2_check == "": # Should see no differences between results and correct files
        print "Full Test Passed"
        full_test.cleanup(other_dirs=[final_dir])
    else:
        print "Trim_Match Test Failed"
        print "sdiff -s test_output/combined_R1.fastq Test_Correct/R1_35bp_test_correct.fq"
        print test_R1_check
        print "sdiff -s test_output/combined_R2.fastq Test_Correct/R2_35bp_test_correct.fq"
        print test_R2_check

if __name__ == '__main__':
    main()
    
