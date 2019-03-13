import os

from Source import CleanGenomic, MapQuantifyGene, Misc, Parameters, Preprocess, QuantifyTail, Results

"""" FLAMSeqMaster():

- Load parameter File
- check/create output dir
- Check/create stats dir

"""

class FLAMSeqMaster():

    def __init__(self, parameterYamlFile):

        print("-----------------")
        print("FLAM-Seq Pipeline")
        print("-----------------")

        # Generate Parameters Object
        self.parameters = Parameters.Parameters(parameterYamlFile)

        # Check outputDir
        expOutDir = self.parameters.getParametersByKey('experiment')['outputDir']
        Misc.Misc.checkDir(expOutDir, 'OutputDir')
        #self.parameters.registerDir('outdir', expOutDir)

        # Check statsDir
        statsOutDir = os.path.join(expOutDir, 'statsDir')
        Misc.Misc.checkDir(statsOutDir, 'statsDir')
        #self.parameters.registerDir('statsdir', statsOutDir)

    def getParameters(self):
        return self.parameters.getParameters()

    def preprocess(self):
    # Preprocess Reads, filter reads with correct adapter and poly(A) tail

        preprocess = Preprocess.Preprocess(self.parameters)
        preprocess.process()

    def quanifyTail(self):
    # Quantify Tail Length from Preprocessed Fastq Files

        quanitfyTail = QuantifyTail.QuantifyTail(self.parameters)
        quanitfyTail.quantifyTails()

    def mapQuantifyGene(self):
    # Map Reads Using StarLong

        mapQuantGene = MapQuantifyGene.MapQuantifyGene(self.parameters)
        mapQuantGene.mapQuantify()

    def genomicCleanTails(self):
    # Remove Genomic Sequence Stretches from polyA tail 5' ends

        cleanGenomic = CleanGenomic.CleanGenomic(self.parameters)
        cleanGenomic.cleanTails()

    def results(self):
    # Aggregate results

        results = Results.Results(self.parameters)
        results.collect()
