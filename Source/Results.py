import os
import pandas as pd
from Source import Misc

"""
    class Results()

    - create geneLengthTable.txt file
    - Aggregate UMI per for each read, polyA Tail Length and sequence and gene info
"""


class Results():

    def __init__(self, parameters):

        print('---------------')
        print('Collect Results')
        print('---------------')

        outputDir = parameters.getParametersByKey('experiment')['outputDir']
        self.sampleName = parameters.getExpName()

        self.resultDir = os.path.join(outputDir, 'resultDir')
        Misc.Misc.checkDir(self.resultDir, 'resultDir')

        # Define input samples
        # - featCnt File
        # - genCleanTail File
        # - umi File

        self.featCntPath = os.path.join(outputDir, 'mapQuantGeneDir', self.sampleName + '_Aligned.sortedByCoord.out.bam.featureCounts')
        self.genCleanTailLengthPath = os.path.join(outputDir, 'cleanGenomicDir', self.sampleName + '_clean_genomic_tail_length.txt')
        self.umiPath = os.path.join(outputDir, 'quantTailDir', self.sampleName + '_umis.txt')

        Misc.Misc.checkFile(self.featCntPath, 'FeatureCountBAM')
        Misc.Misc.checkFile(self.genCleanTailLengthPath, 'Genomic Clean Tail Length')
        Misc.Misc.checkFile(self.umiPath, 'UMI')

    def generateGeneLengthTable(self):
    # Method joins FeatureCount, genomicCleanTail Length and UMI data into one file

        geneLengthTablePath = os.path.join(self.resultDir, self.sampleName + '_gene_polyA_length.csv')

        # Load data as pandas dataFrames
        featCntDf = pd.read_csv(self.featCntPath, sep='\t', header=None)
        featCntDf = featCntDf.iloc[:,[0,3]]
        featCntDf = featCntDf.dropna()
        featCntDf.columns = ['read','gene']

        genCleanTailDf = pd.read_csv(self.genCleanTailLengthPath, sep=',', header=None)
        genCleanTailDf.columns = ['length','tail','read']
        genCleanTailDf['read'] = [rName.lstrip('@') for rName in genCleanTailDf['read']]

        umiDf = pd.read_csv(self.umiPath, sep=',', header=None)
        umiDf.columns = ['umi','read']

        # Merge Dataframes by 'read'

        geneTableMerge = pd.merge(featCntDf, genCleanTailDf, on='read', how='inner')
        geneTableMerge = pd.merge(geneTableMerge, umiDf, on='read', how='inner')
        geneTableMerge = geneTableMerge.loc[:,['read','gene','umi','length','tail']]
        geneTableMerge = geneTableMerge.drop_duplicates(subset='read')
        geneTableMerge.to_csv(geneLengthTablePath, sep=' ', index=False)

    def collect(self):
        self.generateGeneLengthTable()