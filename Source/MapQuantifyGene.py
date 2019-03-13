import os
from Source import Misc
import subprocess


"""
MapReads()

- Map Trimmed Reads using STARLong and generate sorted BAM file


"""


class MapQuantifyGene():

    def __init__(self, parameters):

        print("------------------")
        print("Map Quantify Genes")
        print("------------------")

        outputDir = parameters.getParametersByKey('experiment')['outputDir']
        self.sampleName = parameters.getExpName()

        quantTailDir = os.path.join(outputDir, 'quantTailDir')
        self.sampleTrimFqPath = os.path.join(quantTailDir, self.sampleName + '_trimmed_tail.fq')
        Misc.Misc.checkFile(self.sampleTrimFqPath, 'Fastq Trimmed polyA Tail File')

        # Generate new dir for writing output files
        self.mapQuantGeneDir = os.path.join(outputDir, 'mapQuantGeneDir')
        Misc.Misc.checkDir(self.mapQuantGeneDir, 'mapQuantGeneDir')

        # Check if software is accressible
        self.starLongPath = parameters.getParametersByKey('software')['STARlong']
        self.featureCountPath = parameters.getParametersByKey('software')['featureCounts']
        self.nThreads = parameters.getParametersByKey('software')['nThreads']

        Misc.Misc.checkProg(self.starLongPath, 'STARLong')
        Misc.Misc.checkProg(self.featureCountPath, 'featureCounts')

        # Check STARLong Index
        self.indexDir = parameters.getParametersByKey('experiment')['genomeIndexDir']
        if not os.path.exists(self.indexDir):
            raise FileNotFoundError("Cannot Locate STAR Index at {}\n".format(self.indexDir))

        # Check GTF File
        self.gtfFile = parameters.getParametersByKey('experiment')['annotationGTF']
        Misc.Misc.checkFile(self.gtfFile, 'Annotation GTF')

        # Define Stats dir
        self.statsDir = os.path.join(outputDir, 'statsDir')

    def mapQuantify(self):

        # Map Reads Using StarLong
        self.mapStarLong()
        self.quantifyFeatCounts()
        self.collectStats()

    def mapStarLong(self):

        samplePrefix = os.path.join(self.mapQuantGeneDir, self.sampleName + '_')

        starLongCmd = [self.starLongPath,
                       '--genomeDir', self.indexDir,
                       '--readFilesIn', self.sampleTrimFqPath,
                       '--outFileNamePrefix', samplePrefix,
                       '--runThreadN', str(self.nThreads),
                       '--outFilterMultimapScoreRange', str(20),
                       '--outFilterScoreMinOverLread', str(0),
                       '--outFilterMatchNminOverLread', str(0.66),
                       '--outFilterMismatchNmax', str(1000),
                       '--winAnchorMultimapNmax', str(200),
                       '--seedSearchStartLmax', str(12),
                       '--seedPerReadNmax', str(100000),
                       '--seedPerWindowNmax', str(100),
                       '--alignTranscriptsPerReadNmax', str(100000),
                       '--alignTranscriptsPerWindowNmax', str(10000),
                       '--outSAMtype', 'BAM', 'SortedByCoordinate']

        subprocess.call(starLongCmd)

    def quantifyFeatCounts(self):

        sampleMapBAMPath = os.path.join(self.mapQuantGeneDir, self.sampleName + '_Aligned.sortedByCoord.out.bam')
        sampleFeatCountPath = os.path.join(self.mapQuantGeneDir, self.sampleName + '_fcounts.txt')

        featCountsCmd = [self.featureCountPath,
                         '-L', '-g', 'gene_name', '-s', str(2),
                         '-R', 'CORE',
                         '-a', self.gtfFile,
                         '-o', sampleFeatCountPath,
                         sampleMapBAMPath]

        subprocess.call(featCountsCmd)

    def collectStats(self):
    # Extract relevant stats from Mapping log.final.out and featCount summary
        featCountSumPath = os.path.join(self.mapQuantGeneDir, self.sampleName + '_fcounts.txt.summary')
        starLogPath = os.path.join(self.mapQuantGeneDir, self.sampleName + '_Log.final.out')

        mapQuantStatsPath = os.path.join(self.statsDir, self.sampleName + '_MapQuant.txt')
        mapQuantStats = open(mapQuantStatsPath, 'w')

        with open(featCountSumPath, 'r') as featCountSum, open(starLogPath, 'r') as starLog:
            for line in featCountSum:
                if line.startswith('Unassigned') or line.startswith('Assigned'):
                    mapQuantStats.write('FeatCount_' + line)
            for line in starLog:
                line = line.split('|')
                if not len(line) == 2:
                    continue
                name = line[0].strip()
                num = line[1].strip()

                # Extract Mapping Stats
                if name.startswith('Number of input reads'):
                    mapQuantStats.write('STAR Num Input' + '\t' + num + '\n')
                elif name.startswith('Average input read length'):
                    mapQuantStats.write('STAR Input Read Length' + '\t' + num + '\n')
                elif name.startswith('Uniquely mapped reads number'):
                    mapQuantStats.write('STAR Num Unique Map Reads' + '\t' + num + '\n')
                elif name.startswith('Average mapped length'):
                    mapQuantStats.write('STAR Map Read Lenght' + '\t' + num + '\n')
                elif name.startswith('Number of reads mapped to multiple loci'):
                    mapQuantStats.write('STAR Multimap Reads' + '\t' + num + '\n')

                # Extract Mapping Error Stats
                elif name.startswith('Mismatch rate per base'):
                    num = num.replace('%','')
                    num = float(num) / 100
                    mapQuantStats.write('STAR Mismatch Rate Base' + '\t' + str(num) + '\n')
                elif name.startswith('Deletion rate per base'):
                    num = num.replace('%', '')
                    num = float(num) / 100
                    mapQuantStats.write('STAR Deletion Rate Base' + '\t' + str(num) + '\n')
                elif name.startswith('Insertion rate per base'):
                    num = num.replace('%', '')
                    num = float(num) / 100
                    mapQuantStats.write('STAR Insertion Rate Base' + '\t' + str(num) + '\n')
        mapQuantStats.close()
