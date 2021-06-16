import pysam
import os
from difflib import SequenceMatcher
from Source import Misc

"""
    class CleanGenomic cleans poly(A) tail sequences by comparing 5' ends of poly(A) tail with UTR sequences



"""


class CleanGenomic():

    def __init__(self, parameters):

        print("--------------------------------------------------------")
        print("CleanGenomic: Remove Genomic Sequence from poly(A) Tails")
        print("--------------------------------------------------------")

        # Define Output Dir
        outputDir = parameters.getParametersByKey('experiment')['outputDir']
        self.sampleName = parameters.getExpName()
        self.parameters = parameters

        self.cleanGenDir = os.path.join(outputDir, 'cleanGenomicDir')
        Misc.Misc.checkDir(self.cleanGenDir, 'cleanGenomicDir')

        # Define Input Files
        self.filteredFastqPath = os.path.join(outputDir, 'preprocessDir', self.sampleName + '_preprocessed_filtered.fq')
        self.mappedBAMPath = os.path.join(outputDir, 'mapQuantGeneDir', self.sampleName + '_Aligned.sortedByCoord.out.bam')
        self.genomeFastaPath = parameters.getParametersByKey('experiment')['genomeFasta']
        self.cleanPaLenPath = os.path.join(outputDir, 'quantTailDir', self.sampleName + '_tail_length.txt')

        self.inputFastqPath = parameters.getParametersByKey('experiment')['rawFastq']

        Misc.Misc.checkFile(self.filteredFastqPath, 'FilteredFastq')
        Misc.Misc.checkFile(self.mappedBAMPath, 'MappdBAM')
        Misc.Misc.checkFile(self.genomeFastaPath, 'GenomeFasta')
        Misc.Misc.checkFile(self.cleanPaLenPath, 'CleanPolyA')
        Misc.Misc.checkFile(self.inputFastqPath, 'Input Fastq')

        # Define Parameters for Tail Clearing
        self.alnRawOverlap = 50
        self.genomicSeqExtLength = 300
        self.allowedMismatch = 3
        self.minPaSeqLen = 10

        print('Input Fastq Path: {}'.format(self.filteredFastqPath))
        print('Mapped BAM Path: {}'.format(self.mappedBAMPath))
        print('Genome Fasta Path: {}'.format(self.genomeFastaPath))
        print('PolyA Tail Length Path: {}'.format(self.cleanPaLenPath))

        # Define Stats Dir
        self.statsDir = os.path.join(outputDir, 'statsDir')


    def get_chrom_name(self, tid, bam_aln_file):

        return (bam_aln_file.get_reference_name(tid))

    def matchingString(self, x, y):
        match = ''
        for i in range(0, len(x)):
            for j in range(0, len(y)):
                k = 1
                # now applying while condition untill we find a substring match and length of substring is less than length of x and y
                while (i + k <= len(x) and j + k <= len(y) and x[i:i + k] == y[j:j + k]):
                    if len(match) <= len(x[i:i + k]):
                        match = x[i:i + k]
                    k = k + 1
        return match

    def parseCleanPolyA(self):
    # Parse _tail_length.txt file and return dict with polyA sequence for each read name

        rNamePaSeqDict = dict()

        with open(self.cleanPaLenPath) as paSeqFile:
            for line in paSeqFile:
                _,paSeq,rName = line.rstrip().split(',')
                rNamePaSeqDict[rName] = paSeq
        return rNamePaSeqDict

    # DEPRECATED FOR RELEASE VERSION; causes different results than
    # clean_genomic.py
    # Use for subsequent versions
    def parseFilterFastq(self):
    # Parse _preprocessed_filtered.fq and return dict with rName and Sequence

        rNameFastqSeqDict = dict()

        with pysam.FastxFile(self.filteredFastqPath) as fqFile:
            for e in fqFile:
                rNameFastqSeqDict[e.name] = e.sequence
        return rNameFastqSeqDict


    # Parse Raw Reads for processing
    def parseRawFastq(self):

        rNameInputFqSeqDict = dict()

        # Define RT primer sites in fwd and rev orientation
        adapter = self.parameters.getParametersByKey("experiment")["adapter"]
        ct_fwd_seq = Misc.Misc.rcSeq(adapter)[-4:] + "T" * 4
        ag_rev_seq = "A" * 4 + adapter[:4]

        with pysam.FastxFile(self.inputFastqPath) as inputFq:
            for e in inputFq:
                # Remove Backslash from Fastq Read Names
                e.name = e.name.replace('/', '_')

                eSeq = e.sequence

                if ct_fwd_seq in eSeq[:100]:
                    rNameInputFqSeqDict[e.name] = eSeq
                elif ag_rev_seq in eSeq[-100:]:
                    rNameInputFqSeqDict[e.name] = Misc.Misc.rcSeq(eSeq)

        return rNameInputFqSeqDict

    def cleanTails(self):
    # Method compares 5' ends of poly(A) tail to genomic sequence of respective genes and writes files w/
    # - clean BAM alignments
    # - pa_cleaned_out: poly(A) sequences of each read
    # - non_templated_out: non-poly(A) sequence of each read

        # Get Input data

        rNamePaDict = self.parseCleanPolyA()
        rNameFastqDict = self.parseRawFastq()

        genomeFasta = pysam.FastaFile(self.genomeFastaPath)
        alnBamFile = pysam.AlignmentFile(self.mappedBAMPath, 'rb')

        # Define Output Files
        paGenCleanOutPath = os.path.join(self.cleanGenDir, self.sampleName + '_clean_genomic_tail_length.txt')
        nonTempGenCleanOutPath = os.path.join(self.cleanGenDir, self.sampleName + '_genomic_non_temp_tails.txt')
        cleanedBAMPath = os.path.join(self.cleanGenDir, self.sampleName + '_cleaned.bam')

        # Define out files
        paGenCleanOut= open(paGenCleanOutPath, 'w')
        nonTempGenCleanOut = open(nonTempGenCleanOutPath, 'w')
        cleanedBAM = pysam.AlignmentFile(cleanedBAMPath, 'wb', template = alnBamFile)

        for aln in alnBamFile:

            # Remove non-primary alignements
            if aln.is_unmapped or aln.is_secondary:
                continue

            # Compare query name to raw_seq and pa_seq
            q_name = aln.query_name

            if (q_name in rNameFastqDict) and (q_name in rNamePaDict):

                aln_raw_seq = rNameFastqDict[q_name]
                aln_pa_seq = rNamePaDict[q_name]

                # Get aligned seq part of current read
                query_aln_seq = aln.query_alignment_sequence
                query_aln_seq = query_aln_seq.upper()
                # Define genomic start/end coordinates for read (for + genes: aln_ref_end -> UTR End)
                #                                               (for - genes: aln_ref_start -> UTR End)
                aln_ref_start = aln.reference_start
                aln_ref_end = aln.reference_end

                # Read are default in cDNA orientation and converted to sense orientation by alignment ->
                # RC query alignment seq for + genes, i.e. alignment is reverse

                if aln.is_reverse:
                    query_aln_seq = Misc.Misc.rcSeq(query_aln_seq)

                # Find start coordinate of alignment sequence in raw read and add 50 nt offset
                # First index in raw read which is covered by aligned sequence

                # Raw    [index 0]   5'-Adapter----polyA----UTR.....-3' [index n]
                # Query                                  5'-UTR-...Frac----3'

                aln_raw_read_start = aln_raw_seq.find(query_aln_seq)
                aln_raw_read_start_add = aln_raw_read_start + self.alnRawOverlap

                # If aligned seq is found in raw read: Get genomic sequence and extract templated / non-templated
                # fractions
                if not aln_raw_read_start == -1:

                    # Raw Read Sequence which overlaps with alignment by aln_raw_overlap and contains polyA
                    # and non-mapped fraction of read
                    nmap_raw = aln_raw_seq[:aln_raw_read_start_add]

                    # Get genomic coordiantes for 'non_mapped' sequence.
                    # For + gene nmap_raw proportion starts at alignment end (- aln_raw_overlap) -> aln.is_reverse
                    # For - gene nmap_raw proportion starts at alignment start (+ aln_raw_overlap) ->NOT aln.is_reverse
                    # Retrive genomic sequence with length genomic_seq_ext_length for individual nucleotide mismatch
                    # Mapping

                    if aln.is_reverse:

                        # Get genomic sequence starting 50 nt before end of alignment
                        nmap_genomic_start = aln_ref_end - self.alnRawOverlap
                        try:
                            nmap_genomic_sequence = genomeFasta.fetch(self.get_chrom_name(aln.reference_id, alnBamFile),
                                                                       nmap_genomic_start,
                                                                       nmap_genomic_start + self.genomicSeqExtLength)
                        except KeyError:
                            continue
                        except ValueError:
                            continue
                        nmap_genomic_sequence = nmap_genomic_sequence.upper()
                        nmap_genomic_sequence = Misc.Misc.rcSeq(nmap_genomic_sequence)


                    else:
                        nmap_genomic_start = aln_ref_start + self.alnRawOverlap
                        try:
                            nmap_genomic_sequence = genomeFasta.fetch(self.get_chrom_name(aln.reference_id, alnBamFile),
                                                                       nmap_genomic_start - self.genomicSeqExtLength,
                                                                       nmap_genomic_start)
                        except KeyError:
                            continue
                        except ValueError:
                            continue
                        nmap_genomic_sequence = nmap_genomic_sequence.upper()

                    # Compare nmap_genomic_sequence to nmap_raw sequence and identify position with
                    # 3 mismatches in a row, then jump back and identify last templated nucleotides
                    # Compare nucleotides from 3' end (last index) in reverse order -> This is searching from 3'UTR
                    # towards polyA Tail.

                    # Define last index of raw read 'non-mapped fraction'
                    # Define last index of genomic sequence corr. 'non-mapped fraction'
                    nmap_raw_index_max = len(nmap_raw) - 1
                    nmap_gen_index_max = len(nmap_genomic_sequence) - 1
                    mism_count = 0

                    # Define modulation index for indel cases (i.e. shift seq by n)
                    r_ind_nmap_raw_mod = 0
                    r_ind_nmap_gen_mod = 0

                    # For every nucleotide in nmap_raw (unmapped fraction of raw seq + offset)
                    for ind in range(0, len(nmap_raw)):

                        # Get index and get nucleotide from 3' end on
                        r_ind_nmap_raw = nmap_raw_index_max - ind + r_ind_nmap_raw_mod
                        r_ind_nmap_gen = nmap_gen_index_max - ind + r_ind_nmap_gen_mod

                        # Define nucletide at position for comparison
                        nuc_nmap_raw = nmap_raw[r_ind_nmap_raw]
                        nuc_nmap_genomic_seq = nmap_genomic_sequence[r_ind_nmap_gen]

                        # Check if nucleotides match at position. If, set mism counter 0,
                        # else increase mism counter.

                        if nuc_nmap_raw == nuc_nmap_genomic_seq:
                            mism_count = 0
                        else:
                            mism_count += 1

                        # If three mismatches in a row are reached,
                        # a) Check InDel Case: 1) Check deletion ref gen sequence
                        #                      2) Check insertin ref gen sequence
                        # step back 3 nt and define this as
                        # last templated nucleotide

                        if mism_count == self.allowedMismatch:

                            # Check deletion case. If 3 mismatches are reached, move back 3 ind on r_ind_nmap_raw,
                            # 2 ind on r_ind_nmap_genomic and compare 3 nucleotides at
                            # r_ind_nmap_raw[r_ind_nmap_raw+3:r_ind_nmap_raw] and
                            # r_ind_nmap_gen[r_ind_nmap_gen_ind+2:r_ind_nmap_gen_ind-1]
                            r_ind_nmap_raw_del = r_ind_nmap_raw + self.allowedMismatch
                            r_ind_nmap_gen_del = r_ind_nmap_gen + self.allowedMismatch - 1

                            # Define raw and genomic sequence for deletion case in raw seq at - mism_count postion
                            nmap_raw_del_nuc = nmap_raw[(r_ind_nmap_raw_del - self.allowedMismatch):
                                                        (r_ind_nmap_raw_del + 3 - self.allowedMismatch)]
                            nmap_genomic_del_nuc = nmap_genomic_sequence[(r_ind_nmap_gen_del - self.allowedMismatch):
                                                                         (r_ind_nmap_gen_del + 3 - self.allowedMismatch)]

                            # Check insertion case. If 3 mismatches are reached, move back 2 ind on r_ind_nmap_raw,
                            # 3 ind on r_ind_nmap_genomic and compare 3 nucleotides at
                            # r_ind_nmap_raw[r_ind_nmap_raw-3:r_ind_nmap_raw] and
                            # r_ind_nmap_gen[r_ind_nmap_gen_ind-2:r_ind_nmap_gen_ind+1]
                            r_ind_nmap_raw_ins = r_ind_nmap_raw + self.allowedMismatch - 1
                            r_ind_nmap_gen_ins = r_ind_nmap_gen + self.allowedMismatch

                            # Define raw and genomic sequence for Insertion Case
                            nmap_raw_ins_nuc = nmap_raw[(r_ind_nmap_raw_ins - self.allowedMismatch):
                                                        (r_ind_nmap_raw_ins + 3 - self.allowedMismatch)]
                            nmap_genomic_ins_nuc = nmap_genomic_sequence[(r_ind_nmap_gen_ins - self.allowedMismatch)
                                                                         :(r_ind_nmap_gen_ins + 3 - self.allowedMismatch)]

                            # If deletion is defined, set modulation_index for raw_seq and genomic_seq
                            if nmap_raw_del_nuc == nmap_genomic_del_nuc:
                                # Set modulation indices and set mism_count = 0
                                r_ind_nmap_raw_mod += self.allowedMismatch
                                r_ind_nmap_gen_mod += self.allowedMismatch - 1
                                mism_count = 0
                            # If insertion is defined, set modulation_index for raw_seq and genomic_seq
                            elif nmap_raw_ins_nuc == nmap_genomic_ins_nuc:
                                # Set modulation indices and set mism_count = 0
                                r_ind_nmap_raw_mod += self.allowedMismatch - 1
                                r_ind_nmap_gen_mod += self.allowedMismatch
                                mism_count = 0
                            # If neither deletion nor insertion detected, define indices for genomic seq and raw seq templated
                            # and non-templated
                            else:
                                last_index_nmap_raw = r_ind_nmap_raw + self.allowedMismatch
                                last_index_nmap_genomic = r_ind_nmap_gen + self.allowedMismatch
                                break
                    else:
                        # If loop reaches end, continue with next alignment
                        continue

                    # Define raw_non_templated fraction of raw reads, these should contain polyA Tail and adaptors
                    nmap_raw_non_templated = nmap_raw[:last_index_nmap_raw]

                    # Compare 5' Trimmed Tail (nmap_raw_non_templated) to orginial Tail sequence
                    # Find longest common substring between strings
                    # ! Use autojunk=False parameter to prevent uncontrolled substring identification for long substrings

                    pa_substr_matcher = SequenceMatcher(None, nmap_raw_non_templated, aln_pa_seq, autojunk=False)
                    pa_substr_match = pa_substr_matcher.find_longest_match(0, len(nmap_raw_non_templated), 0,
                                                                           len(aln_pa_seq))

                    pa_substr = nmap_raw_non_templated[pa_substr_match.a:(pa_substr_match.a + pa_substr_match.size)]

                    # If len of pa_substring >= min_pa_seq_len, write pA Tail, len and r_name. Additionally write 'correct' alignments to
                    # new BAM file

                    if len(pa_substr) >= self.minPaSeqLen:
                        # Write genomic_clean-pa_len to file
                        pa_cor_line = ','.join([str(len(pa_substr)), pa_substr, '@' + q_name]) + '\n'
                        paGenCleanOut.write(pa_cor_line)
                        cleanedBAM.write(aln)

                        # Write non_templated_read to file
                        non_temp_cor_line = ','.join(
                            [str(len(nmap_raw_non_templated)), nmap_raw_non_templated, '@' + q_name]) + '\n'
                        nonTempGenCleanOut.write(non_temp_cor_line)

        paGenCleanOut.close()
        nonTempGenCleanOut.close()
        alnBamFile.close()
