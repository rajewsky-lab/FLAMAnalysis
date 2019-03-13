import argparse
import FLAMSeqMaster

#########
# about #
#########

# __version__ = "1.0.0"
# __author__ = ["Nikos Karaiskos", "Jonathan Alles", "Ivano Legnini", "Salah Ayoub"]
# __status__ = "beta"
# __licence__ = "GPLv3"
# __email__ = {"nikolaos.karaiskos", "jonathan.alles", "ivano.legnini", "salah.ayoub"}@mdc-berlin.de



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Process FLAM-Seq Sequencing Data')
    subparser = parser.add_subparsers(dest="command")
    parser_all = subparser.add_parser('all', help="Run Complete FLAM-Seq Pipeline")
    parser_all.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    parser_preprocess = subparser.add_parser('preprocess', help="Preprocess FLAM-Seq Fastq")
    parser_preprocess.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    parser_quantTail = subparser.add_parser('quantTail', help="Quantify Tail Length for Preprocessed Fastq")
    parser_quantTail.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    parser_mapQuant = subparser.add_parser('mapQuant', help="Map Reads using STARLong and Quantify Genes")
    parser_mapQuant.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    parser_cleanGenomic = subparser.add_parser('cleanGenomic',
                                               help="Compare poly(A) Tails to Genomic Sequence to Clean Tail Sequence")
    parser_cleanGenomic.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    parser_results = subparser.add_parser('result',help="Collect polyA tail length estimates and expression data")
    parser_results.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

    params = parser.parse_args()

    command = params.command
    parameters = params.parameters

    if command == 'all':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.preprocess()
        master.quanifyTail()
        master.mapQuantifyGene()
        master.genomicCleanTails()
        master.results()

    elif command == 'preprocess':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.preprocess()

    elif command == 'quantTail':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.quanifyTail()

    elif command == 'mapQuant':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.mapQuantifyGene()

    elif command == 'cleanGenomic':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.genomicCleanTails()

    elif command == 'result':
        master = FLAMSeqMaster.FLAMSeqMaster(parameterYamlFile=parameters)

        master.results()

    else:
        print('Specify Correct Command')