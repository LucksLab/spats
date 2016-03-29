#!/usr/bin/python

#Brief module for reverse-complementing a fastq file with fastx_reverse_complement
#Designed to be replaced with code that does fast revcomping to remove dependency on fastx tools

import os


def rev_comp(inputfile, output_revcomp='blank'):
    
    if output_revcomp == 'blank':
        output_revcomp = os.path.splitext(inputfile)[0] + '_rc.fq'
    
   #Uses the given input and output filenames to make a reverse complement, outputs file name
    os.system("fastx_reverse_complement -Q33 -i {0} -o {1}".format(inputfile,output_revcomp))


def main(argv=None,):
   pass 
    


if __name__ == "__main__":
    main()
    sys.exit(0)