# FLAM-Seq Data Analysis Pipeline

Welcome to FLAM-Seq, nice to see you again.

## Software Dependencies

The analysis pipeline is implemented in Python3 and requires <a href="https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64/STARlong">STARlong</a>
and <a href="http://subread.sourceforge.net/">featureCounts</a> as addtional software. Packages and Software Versions used are listed below:

```
python 3.6.7 Anaconda, Inc.
STARlong 2.5.4b
FeatureCounts 1.6.0
regex 2018.2.21
pysam 0.14
pandas 0.23.4
yaml 0.1.7
```

## Usage


    -h:                 Show help
    -p, --parameters:   Parameter YAML file specifying configurations for analysis

    all                 Run complete FLAM-Seq Analysis Pipeline. This includes all steps listed below:
                        preprocess, quantTail, mapQuant, cleanGenomic, result

    preprocess          Preprocess input fastq reads. Check if reads contain adapters and place reads in correct orientation.
                        This step creates /preprocessDir and creates *_preprocessed_filtered.fq files
                        containing preprocessed reads and *_preprocessed_filtered_error.fq, containing reads without
                        adapters/poly(A) tails.

    quantTail           Quantify poly(A) tail length and sequence from preprocessed reads. This step uses 2 Algorithms
                        with different parameter combinations in order to extract putative poly(A) tail sequences from
                        reads. This step creates /quantTailDir with *_trimmed_tail.fq, containing reads with
                        poly(A) tail trimmed for later mapping, *_umis.txt, containing read names and read UMIs,
                        *_tail_length.txt containing poly(A) tail length, sequence and read name and
                        *_no_tail.fq containing reads without detectable poly(A) sequence.

    mapQuant            Map trimmed reads and assign reads to genes. This step maps reads with trimmed poly(A) tail
                        sequences using STARLong to genome index specified in parameters.yaml. Next aligned reads
                        are assigned to gene models using FeatureCounts and GTF file specified in parameters.yaml.
                        This step creates /mapQuantGeneDir with *__Aligned.sortedByCoord.out.bam containing
                        mapped read BAM file (sorted), *_Aligned.sortedByCoord.out.bam.featureCounts containing read names
                        tagged with gene-of-origin.

    cleanGenomic        cleanGenomic compares the sequences derived from putative poly(A) tail with the genomic sequence
                        of the mapping location of the respective read. This step 'cleans' nucleotides from 5'ends of putative
                        poly(A) sequences from nucleotides which are encoded in the genome.
                        This step requires genome fasta file specified in parameters.yaml.
                        It creates /cleanGenomicDir with *_cleaned.bam conataining aligned reads filtered by
                        reads with valid poly(A) tail, *_clean_genomic_tail_length.txt, containing cleaned poly(A)
                        length/sequence and *_genomic_non_temp_tails.txt containing parts of each read that are not
                        templated by genome.

    result              Aggregate data from above steps into *_gene_polyA_length.csv. This file contains for each
                        read poly(A) tail length and sequence, UMI and gene the read maps to. This step creates
                        /resultDir.


## Perform FLAM-Seq Analysis

Before running this analysis pipeline, PacBio or Nanopore sequencing data need to ne converted to ```.fastq``` format.
This can be done using default software of respective sequencing systems.
The obtained ```.fastq``` file is the starting point for the analysis.

Running the pipeline requires only configuration of the ```.parameters.yaml``` file. This file specifies all relevant
parameters for the analysis, such as path to input files, output directory, mapping index, etc.

It is recommended to generate a new ```parameters.yaml``` file for each sample to be analyzed and store it in the
respective output folder.

The pipeline can the simply be run using the ```FLAMSeqAnalysis.py``` script:

```python3 /path/to/FLAMSeqAnalysis.py command -p /path/to/parameter.yaml```


e.g. for running the complete pipeline

```python3 /path/to/FLAMSeqAnalysis.py all -p /path/to/parameter.yaml```

or for mapping and quantifying the reads

```python3 /path/to/FLAMSeqAnalysis.py mapQuant -p /path/to/parameter.yaml```



## Parameters.yaml

A template ```parameters.yaml``` file is provided with the pipeline. The file can be renamed but structure of keywords withing the
file (everything before the ':') must not be changed. All paths need to be adapted to user requirements. An example is
shown below:


```
experimentName: "FLAMSeq_1"               # Name of Experiment. This will be prefix for generated analysis files.

experiment:
  rawFastq: "/path/to/pacbioreads.fq"     # Define Path to Input PacBio / Nanopore Reads in fastq format
  outputDir: "/path/to/outputDir"         # Define Path to Output dir for writing analysis results. Pipeline will generate of dir is not present.
  genomeIndexDir: "/path/to/STARIndex"    # Define Dir for STARLong Index for Mapping
  annotationGTF: "/path/to/GTFFile"       # GTF File containing gene annotations for assigning reads to genes
  genomeFasta: "/path/to/genomeFasta"     # Genome Fasta File containing genome sequence (that the index was generated from)

software:
  STARlong: "/path/to/STARlong"           # Path to STARlong executable
  featureCounts: "/path/to/featureCounts" # Path to FeatureCounts executable
  nThreads: 8                             # Number of threads for mapping reads
```

## Results

We recommend continuing your analysis with the ```*_gene_polyA_length.csv``` file generated in the ```/resultDir```
after running the complete pipeline. This file can easily be loaded into R or Pandas for downstream analysis on poly(A)
tails.


### Further Info

FLAM-Seq has beed developed by Ivano Legnini, Salah Ayoub, Nikos Karaiskos and Jonathan Alles as members of the
<a href= "https://www.mdc-berlin.de/n-rajewsky">Nikolaus Rajewsky Lab</a> at the <a href="https://www.mdc-berlin.de/">
Max Delbruck Center Berlin</a>. FLAM-Seq allows for sequencing of entire RNA molecules. A preprint describing the FLAM-Seq
method can be found <a href="https://www.biorxiv.org/content/10.1101/547034v1">here</a>.
