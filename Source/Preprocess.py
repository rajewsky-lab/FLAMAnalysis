import pysam
from Source import Misc
import os
from itertools import islice

import regex
import sys

"""

class Preprocess()

- generate preprocess dir

- fastqFile w/ oriented reads containing adapters
- fastqFile containing reads w/o adapters/erroneous reads
- stats file containing preprocessing stats
    
"""


# Declare Module Level Vars for defining threshold for defining poly(A) Tails
thresholdA = 10
#thresholdG = 9

#combAG = 'A' * thresholdA + 'G' * thresholdG
#combCT = 'C' * thresholdG + 'T' * thresholdA

class Preprocess():

    # Define Parameters for

    def __init__(self, parameters):
        self.fastqPath = parameters.getParametersByKey('experiment')['rawFastq']
        self.parameters = parameters
        self.sampleName = self.parameters.getExpName()

        # Check if outdir is registered, if so check preprocessDir
        expOutdir = self.parameters.getParametersByKey('experiment')['outputDir']

        # Generate preprocessing dir
        self.preprocessDir = os.path.join(expOutdir, 'preprocessDir')
        Misc.Misc.checkDir(self.preprocessDir, 'preprocessDir')

        # Check input Fastq File
        Misc.Misc.checkFile(self.fastqPath, 'Input Fastq File')

        # Define Stats dir
        self.statsDir = os.path.join(expOutdir, 'statsDir')

    def process(self):

        print('----------------')
        print("Preprocess Reads")
        print("----------------")

        # Create fq files with filtered reads
        prep_dir = self.preprocessDir

        fq_filter_pos_path = os.path.join(prep_dir, self.sampleName + '_preprocessed_filtered.fq')
        fq_filter_error_path = os.path.join(prep_dir, self.sampleName + '_preprocessed_filtered_error.fq')

        fq_filter_pos = open(fq_filter_pos_path, 'w')
        fq_filter_error = open(fq_filter_error_path, 'w')

        # Define Counter for each class of reads
        cnt_filter = 0
        cnt_filter_error = 0

        # Load fastq as pysam alignment
        with pysam.FastxFile(self.fastqPath) as fq_file:

            for cnt, e in enumerate(islice(fq_file,None)):

                if (cnt % 100000 == 0):
                    sys.stdout.write("Processed {} reads\n".format(cnt))

                # Reformat read names  and remove /
                e.name = e.name.replace('/','_')

                # Check if sequence contains correct GI Tail and oligo dT sequence
                res = self.assessRead(e.sequence)

                if res == 1:
                    e.sequence = Misc.Misc.rcSeq(e.sequence)
                    fq_filter_pos.write(str(e) + '\n')
                    cnt_filter = 1
                elif res == -1:
                    fq_filter_pos.write(str(e) + '\n')
                    cnt_filter += 1
                else:
                    fq_filter_error.write(str(e) + '\n')
                    cnt_filter_error += 1

        fq_filter_pos.close()
        fq_filter_error.close()

        cnt_total = cnt_filter + cnt_filter_error

        # Write Stats
        preprocessStatsPath = os.path.join(self.statsDir, self.sampleName + '_Preprocess.txt')
        preprocessStats = open(preprocessStatsPath, 'w')
        preprocessStats.write("Preprocessed Reads Total\t{}\n".format(cnt_total))
        preprocessStats.write("Preprocessed Reads Filtered\t{}\n".format(cnt_filter))
        preprocessStats.write("Preprocessed Reads Error\t{}\n".format(cnt_filter_error))
        preprocessStats.close()


    def assessRead(self, read):
    # Check if reads contains adapter sequence (9xG's) and 10 nt oligo dT or rev comp
        adapter = self.parameters.getParametersByKey("experiment")["adapter"]
        combAG = "A" * thresholdA + adapter
        combCT = Misc.Misc.rcSeq(adapter) + "T" * thresholdA

        if (regex.findall('(' + combAG + '){s<=1}', read) or
                regex.findall('(' + 'AA' + combAG + '){s<=2}', read) or
                regex.findall('(' + 'AAAAA' + combAG + '){s<=3}', read) or
                regex.findall('(' + 'AAAAAAAAAAAAA' + combAG + '){s<=4}', read)):
            return 1
        if (regex.findall('(' + combCT + '){s<=1}', read) or
                regex.findall('(' + combCT + 'TT' + '){s<=2}', read) or
                regex.findall('(' + combCT + 'TTTTT' + '){s<=3}', read) or
                regex.findall('(' + combCT + 'TTTTTTTTTTTTT' + '){s<=4}', read)):
            return -1
        else:
            return 0




