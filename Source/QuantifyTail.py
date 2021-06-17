import os
from Source import Misc

import regex
from math import ceil
from collections import Counter
import pysam
from itertools import islice

#pcr_handle = 'GGTAATACGACTCACTATAGCGAGA'
#pcr_handle2 = 'TGAGTCGGCAGAGAACTGGCGAA'

class QuantifyTail():

    def __init__(self, parameters):
        outputDir = parameters.getParametersByKey('experiment')['outputDir']
        self.sampleName = parameters.getExpName()
        self.parameters = parameters

        # Define Input Reads
        preprocessDir = os.path.join(outputDir, 'preprocessDir')
        self.preprocessedReads = os.path.join(preprocessDir, self.sampleName + '_preprocessed_filtered.fq')
        Misc.Misc.checkFile(self.preprocessedReads, 'preprocessed Fastq file')

        # Define OutFile
        self.quantTailDir = os.path.join(outputDir, 'quantTailDir')
        Misc.Misc.checkDir(self.quantTailDir, 'QuantTailDir')

        # Define StatsDir
        # Define Stats dir
        self.statsDir = os.path.join(outputDir, 'statsDir')

    def quantifyTails(self):

        print("\nQuantify Tail Length")
        print("Input Fastq: {}".format(self.preprocessedReads))

        # Define output files for sorting out reads
        no_tail_fq_path = os.path.join(self.quantTailDir, self.sampleName + '_no_tail.fq')
        trimmed_fq_path = os.path.join(self.quantTailDir, self.sampleName + '_trimmed_tail.fq')
        umi_path = os.path.join(self.quantTailDir, self.sampleName + '_umis.txt')
        tail_length_path = os.path.join(self.quantTailDir, self.sampleName + '_tail_length.txt')

        no_tail_fq = open(no_tail_fq_path, 'w')
        trimmed_fq = open(trimmed_fq_path, 'w')
        umi_txt = open(umi_path, 'w')
        tail_length_txt = open(tail_length_path, 'w')

        cntPaFound = 0
        cntNoPa = 0

        with pysam.FastxFile(self.preprocessedReads) as fq_in:

            for cnt, e in enumerate(islice(fq_in,None)):

                cnt += 1
                if cnt % 10000 == 0:
                    print("Quantified Tails for {} Reads".format(cnt))

                # Check if poly(A) Tail can be found
                # find_tail_length returns tail length, tail sequence, UMI sequence

                tail_gen = self.find_tail_length(e.sequence, threshold=25)

                # If a) tails is found and b) more than 20 nt sequence are left (for later mapping)
                if len(tail_gen[1]) > 0 and len(e.sequence) - (e.sequence.index(tail_gen[1]) + len(tail_gen[1])) >= 20:

                    # Call majority vote func: Run 2 Algorithms for pA Lenght estimation with diff params for quantification
                    res = self.majority_vote(e.sequence)

                    # Sort Reads as:
                    #   a) _trimmed.fastq
                    # -> If no tail is detected
                    #   b) _no_tail.fastq
                    #   c) _umis.txt
                    #   d) _tail_length.txt

                    if (res[0] == 0):
                        no_tail_fq.write(str(e) + '\n')
                        cntNoPa += 1
                    # If poly(A) Tail was found for read:
                    # - Write Fastq w/o tail sequence for mapping
                    # - Write poly(A) length and sequence
                    # - Write UMI and read names

                    else:
                        pa_len = res[0]
                        pa_seq = res[1]
                        umi = res[2]

                    # Remove poly(A) Tail and adapters from sequence and quality strings
                    # Define Last index of adapter + poly(A) sequence
                        pa_end_index = e.sequence.index(pa_seq) + pa_len
                        e.sequence = e.sequence[pa_end_index:]
                        e.quality = e.quality[pa_end_index:]
                        trimmed_fq.write(str(e) + '\n')

                        umi_str = umi + ',' + e.name + '\n'
                        umi_txt.write(umi_str)

                        pa_len_str = ','.join([str(pa_len), pa_seq, e.name]) + '\n'
                        tail_length_txt.write(pa_len_str)
                        cntPaFound += 1

                else:
                    no_tail_fq.write(str(e) + '\n')
                    cntNoPa += 1

        no_tail_fq.close()
        trimmed_fq.close()
        umi_txt.close()
        tail_length_txt.close()

        # Write Stats
        quantTailStatsPath = os.path.join(self.statsDir, self.sampleName + '_QuantifyTails.txt')
        quantTailStats = open(quantTailStatsPath, 'w')
        quantTailStats.write('Total Reads\t{}\n'.format(cntPaFound+cntNoPa))
        quantTailStats.write('Reads with polyA tail\t{}\n'.format(cntPaFound))
        quantTailStats.write('Reads no polyA tail\t{}\n'.format(cntNoPa))
        quantTailStats.close()


    def find_tail_length(self, read, threshold=30):
    # Define poly(T) tail length for sequence using T extension with increased mismatches. Extract UMI sequence.
    # Return tail length, tail sequence and UMI sequence

        # Define init params for T sequence search
        thr = '1'
        polyA = 'T' * 10
        tail = ''

        if regex.findall('(' + polyA + '){s<=1}', read):
            start_pos = read.find(regex.findall('(' + polyA + '){s<=1}', read)[0])
            polyA_found = True
            polyA_length = len(polyA)
            while polyA_found:
                thr = str(len(polyA) // threshold + 1)
                polyA += 'T'
                if regex.findall('(' + polyA + '){s<=' + thr + '}', read[start_pos:]):
                    polyA_length += 1
                else:
                    polyA_found = False
            tail = regex.findall('(' + polyA[:-1] + '){s<=' + thr + '}', read[start_pos:])[0]
            tail = tail[tail.find('T'):]

        if tail == '':
            return 0, '', ''

        # Trim Read until only T is left

        for i in range((len(tail) // threshold) + 1):
            tail = self.chop_tail(tail, -ceil(threshold / 10))

        tail_start = read.find(tail)

        adapter = Misc.Misc.rcSeq(self.parameters.getParametersByKey("experiment")["adapter"])
        pcr_handle = Misc.Misc.rcSeq(self.parameters.getParametersByKey("experiment")["primer3"])

        # UMI extraction
        if read[:tail_start].endswith(adapter):
            umi_end = tail_start - len(adapter)
            umi = read[(umi_end - 10):umi_end]
        elif adapter[:-1] in read[tail_start - len(adapter):tail_start]:
            umi_end = tail_start - len(adapter) + read[tail_start - len(adapter):tail_start].find(adapter[:-1])
            umi = read[(umi_end - 10):umi_end]
        elif read.find(pcr_handle[-10:]) != -1:
            umi = read[read.find(pcr_handle[-10:]) + 10:read.find(pcr_handle[-10:]) + 20]
        else:
            tail = ''
            umi = ''

        return len(tail), tail, umi

    def chop_tail(self, tail, r):
        """Chop a tail on the right end if it shows bases other than T
        tail -- the tail as a string
        r    -- the window, has to be negative"""
        if r > len(tail): # catch exception
            return tail
        first_pos = 0
        for idx in range(-1, r, -1):
            try:
                if tail[idx] != 'T':
                    first_pos = idx
            except IndexError:
                return tail
        if first_pos == 0:
            return tail
        return tail[:first_pos]

    def pa_slide(self, seq, windowsize=20, threshold=0.8):
        t_seq = None
        len_t_seq = 0
        for ind in range(0, (len(seq) - windowsize)):
            substr = seq[ind:(ind + windowsize)]
            # Define T content of substring
            t_cnt = substr.count('T')
            t_perc = t_cnt / windowsize
            if t_perc < threshold:
                # Define leftmost TT dinucleotide
                t_end = substr.rfind('TT')
                if not t_end is None:
                    t_end += ind + 2
                    t_seq = seq[0:t_end]
                    break
                else:
                    continue

        # Chop the last bases if they are not enriched in Ts
        if t_seq != None:
            for i in range((len(t_seq) // 30) + 1):
                t_seq = self.chop_tail(t_seq, -ceil(4))
            len_t_seq = len(t_seq)

        return len_t_seq, t_seq


    def majority_vote(self, read):
        res_dict = {0: self.find_tail_length(read, threshold=25),
                    1: self.find_tail_length(read, threshold=30),
                    2: self.find_tail_length(read, threshold=35),
                    3: self.find_tail_length(read, threshold=40)}

        int1 = read.find(res_dict[0][1])
        int1 += read[int1:].find('T')

        # Change beginning of read
        read = read[int1:]

        pa_dict = {4: self.pa_slide(read, windowsize=20, threshold=0.8),
                   5: self.pa_slide(read, windowsize=20, threshold=0.85),
                   6: self.pa_slide(read, windowsize=25, threshold=0.8),
                   7: self.pa_slide(read, windowsize=25, threshold=0.85),
                   8: self.pa_slide(read, windowsize=30, threshold=0.8),
                   9: self.pa_slide(read, windowsize=30, threshold=0.85)}

        tail_lengths = [res_dict[0][0], res_dict[1][0], res_dict[2][0], res_dict[3][0],
                        pa_dict[4][0], pa_dict[5][0], pa_dict[6][0], pa_dict[7][0],
                        pa_dict[8][0], pa_dict[9][0]]
        c = Counter(tail_lengths)
        value, count = c.most_common()[0]

        idx_value = tail_lengths.index(value)
        if idx_value <= 3:
            return res_dict[idx_value]
        else:
            # This need to be corrected, umi might not be correct
            return pa_dict[idx_value][0], pa_dict[idx_value][1], res_dict[0][2]
