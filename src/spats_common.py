
import os
import string

rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

class FastqRecord(object):

    def __init__(self):
        self.recordNumber = 0
        self.reset()

    def reset(self):
        self.identifier = None
        self.sequence = None
        self.identifier2 = None
        self.quality = None

    def read(self, infile):
        first = infile.readline()
        if not first or len(first) == 0:
            self.reset()
            return
        self.identifier = first.rstrip('\r\n')
        self.sequence = infile.readline().rstrip('\r\n')
        self.identifier2 = infile.readline().rstrip('\r\n')
        self.quality = infile.readline().rstrip('\r\n')

    def write(self, outfile):
        for line in [ self.identifier, self.sequence, self.identifier2, self.quality ]:
            outfile.write(line)
            outfile.write('\n')

    def reverse_complement(self):
        self.sequence = self.sequence.translate(rev_comp_complementor)[::-1]
        self.quality = self.quality[::-1]

def fastq_transform(input_path, output_path, record_transformer):
    with open(input_path, 'rb') as infile:
        with open(output_path, 'wb') as outfile:
            record = FastqRecord()
            while True:
                record.read(infile)
                if not record.identifier:
                    break
                record_transformer(record)
                if record.identifier:
                    record.write(outfile)

def rev_comp(inputfile, output_revcomp = None):
    output_path = output_revcomp or os.path.splitext(inputfile)[0] + '_rc.fq'
    fastq_transform(inputfile, output_path, lambda x : x.reverse_complement())
