
import math
import os
import string


rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

nuc_A = (1)
nuc_C = (1<<1)
nuc_G = (1<<2)
nuc_T = (1<<3)

char_to_mask = { 'A' : nuc_A,
                 'C' : nuc_C,
                 'G' : nuc_G,
                 'T' : nuc_T,
                 'U' : nuc_T,
                 'Y' : (nuc_C | nuc_T),
                 'R' : (nuc_G | nuc_A),
                 'S' : (nuc_C | nuc_G),
                 'W' : (nuc_A | nuc_T),
                 'K' : (nuc_G | nuc_T),
                 'M' : (nuc_A | nuc_C),
                 'B' : (nuc_C | nuc_G | nuc_T),
                 'D' : (nuc_A | nuc_G | nuc_T),
                 'H' : (nuc_A | nuc_C | nuc_T),
                 'V' : (nuc_A | nuc_C | nuc_G),
                 'N' : (nuc_A | nuc_C | nuc_G | nuc_T),
}

def matches_mask(seq, maskvals):
    # assume len(seq) >= len(maskvals) -- otherwise this will raise
    for i in range(len(maskvals)):
        seqval = char_to_mask[seq[i]]
        if 0 == (seqval & maskvals[i]):
            return False
    return True

class Mask(object):

    def __init__(self, chars):
        self.chars = chars
        self.values = [ char_to_mask[ch] for ch in chars ]


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

    def relabel_reads(self):
        self.recordNumber += 1
        self.identifier2 += self.identifier[1:]
        self.identifier = "@{}".format(self.recordNumber)

    def strip_and_relabel(self, size, total):
        self.identifier = "@{}".format(total)
        self.sequence = self.sequence[size:] if size > 0 else self.sequence[:size]
        self.quality = self.quality[size:] if size > 0 else self.quality[:size]
        self.identifier2 = "+"


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

def relabel_reads(handle_reads_path, nonhandle_reads_path, output_path):
    transformer = lambda x : x.relabel_reads()
    fastq_transform(handle_reads_path, os.path.join(output_path, "NOMASK_1.fq"), transformer)
    fastq_transform(nonhandle_reads_path, os.path.join(output_path, "NOMASK_2.fq"), transformer)

def filter_reads(handle_reads_path, nonhandle_reads_path, mask_strings, output_folder, relabel = True):
    total = 0
    num_chucked = 0
    masks = [ Mask(m) for m in mask_strings ]
    maskSize = max([ len(m.chars) for m in masks ])
    with open(handle_reads_path, 'rb') as handle_in:
        with open(nonhandle_reads_path, 'rb') as nonhandle_in:
            for mask in masks:
                mask.handle_out = open(os.path.join(output_folder, "{}_1.fq".format(mask.chars)), 'wb')
                mask.nonhandle_out = open(os.path.join(output_folder, "{}_2.fq".format(mask.chars)), 'wb')
            try:
                handle_record = FastqRecord()
                nonhandle_record = FastqRecord()

                while True:
                    handle_record.read(handle_in)
                    nonhandle_record.read(nonhandle_in)
                    if not handle_record.identifier or not nonhandle_record.identifier:
                        break
                    total += 1
                    matched_mask = None
                    if min(len(handle_record.sequence), len(nonhandle_record.sequence)) > maskSize + 1:
                        for mask in masks:
                            if matches_mask(handle_record.sequence, mask.values):
                                matched_mask = mask
                                break
                    if not matched_mask:
                        num_chucked += 1
                    else:
                        sz = len(mask.chars)
                        handle_record.strip_and_relabel(sz, total)
                        handle_record.write(matched_mask.handle_out)
                        nonhandle_record.strip_and_relabel(-sz, total)
                        nonhandle_record.write(matched_mask.nonhandle_out)
            finally:
                for mask in masks:
                    mask.handle_out.close()
                    mask.nonhandle_out.close()

    print "Kept {} of {} reads".format(total - num_chucked, total)

class SamRecord(object):

    def parse(self, line):
        bits = line.split("\t")
        if len(bits) < 6:
            self.identifier = None
            return
        self.identifier = bits[0]
        self.flags = int(bits[1])
        # bits[2] is the target name
        self.left = int(bits[3]) - 1 # TODO: subtract 1 is required, explain
        self.quality = int(bits[4])
        lengthPart = bits[5][:-1]
        self.length = int(lengthPart if len(lengthPart) > 0 else 0)
        if self.length > 0 and "M" != bits[5][-1]:
            raise Exception("Expected M on SAM length field, got: {}".format(bits[5]))
        # rest of bits are not used afaict
        self.right = self.left + self.length

def fasta_parse(target_path):
    names = []
    seqs = []
    with open(target_path, 'rb') as infile:
        while True:
            line = infile.readline()
            if 0 == len(line):
                break
            names.append(line.strip('>\n'))
            line = infile.readline()
            seqs.append(line.strip())
    return names, seqs

def get_counts(sam_path, target_length):
    num_sites = target_length + 1
    counts = [ 0 for x in range(num_sites) ] # TODO: numpy.empty(num_sites, dtype=int) is better
    num_fragments = 0
    num_chucked = 0
    r1 = SamRecord()
    r2 = SamRecord()
    with open(sam_path, 'rb') as infile:
        while True:
            line = infile.readline()
            if 0 == len(line):
                break
            if line.startswith('@'):
                continue
            r1.parse(line)
            r2.parse(infile.readline())
            if not r1.identifier or r1.identifier != r2.identifier:
                raise Exception("Parse error?") # might just want continue?
            num_fragments += 1
            left = min(r1.left, r2.left)
            right = max(r1.right, r2.right)
            if left >= 0 and left < num_sites and right == target_length:
                # JBL - only register this fragment if the left read is within the sequence and the right read aligns with the end of the RNA (or RNA subsequence)
                counts[left] += 1
            else:
                num_chucked += 1
    return counts, num_fragments, num_fragments - num_chucked

def compute_from_counts(num_sites, treated_counts, untreated_counts):
    n = int(num_sites)
    betas = [ 0 for x in range(n+1) ]
    thetas = [ 0 for x in range(n+1) ]
    treated_sum = 0.0    # keeping a running sum for both
    untreated_sum = 0.0  # channels is much faster
    running_c_sum = 0.0  # can also do it for c

    for k in range(n):
        X_k = float(treated_counts[k])
        Y_k = float(untreated_counts[k])
        treated_sum += X_k    #treated_sum = float(sum(treated_counts[:(k + 1)]))
        untreated_sum += Y_k  #untreated_sum = float(sum(untreated_counts[:(k + 1)]))
        if 0 == treated_sum  or  0 == untreated_sum:
            betas[k] = 0
            thetas[k] = 0
        else:
            Xbit = (X_k / treated_sum)
            Ybit = (Y_k / untreated_sum)
            if Ybit >= 1:
                betas[k] = 0
                thetas[k] = 0
            else:
                betas[k] = max(0, (Xbit - Ybit) / (1 - Ybit))
                thetas[k] = math.log(1.0 - Ybit) - math.log(1.0 - Xbit)
                running_c_sum -= math.log(1.0 - betas[k])

    c = running_c_sum
    c_factor = 1.0 / c
    for k in range(n+1):
        thetas[k] = c_factor * thetas[k]

    return betas, thetas, c

def compute_profiles(target_path, treated_path, untreated_path, output_dir):
    names, seqs = fasta_parse(target_path)

    # TODO: handle multiple sequences?

    name = names[0]
    seq = seqs[0]
    n = len(seq)
    treated_counts, treated_total, treated_kept = get_counts(treated_path, n)
    untreated_counts, untreated_total, untreated_kept = get_counts(untreated_path, n)
    print "Processed {} properly paired fragments, kept {}/{} ({:.2f}%) treated, {}/{} ({:2f}%) untreated".format(treated_total + untreated_total,
                                                                                                                  treated_kept,
                                                                                                                  treated_total,
                                                                                                                  (100.0 * float(treated_kept)) / float(treated_total),
                                                                                                                  untreated_kept,
                                                                                                                  untreated_total,
                                                                                                                  (100.0 * float(untreated_kept)) / float(untreated_total))

    betas, thetas, c = compute_from_counts(n, treated_counts, untreated_counts)

    out_path = os.path.join(output_dir, "reactivities.out")
    with open(out_path, 'wb') as outfile:
        outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
        format_str = "{name}\t{rt}\t".format(name = name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = c)
        # TODO: xref https://trello.com/c/OtbxyiYt/23-3-nt-missing-from-reactivities-out
        # looks like we may want this to be range(n), chopping off was unintentional bug of previous version
        for i in range(n - 1):
            outfile.write(format_str.format(i = i,
                                            nuc = seq[i - 1] if i > 0 else '*',
                                            tm = treated_counts[i],
                                            um = untreated_counts[i],
                                            b = betas[i] if i > 0 else '-',
                                            th = thetas[i] if i > 0 else '-'))
