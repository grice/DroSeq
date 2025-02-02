# Droseq mirbase pipeline

REQUIREMENTS python3

Command line programs:
 - sickle
 - bowtie2
 - samtools
 - FLASH

Python modules:
 - pysam

Bowtie2 database containing pri-miRNA sequences in reverse complement
based on the library construction strategy. This step only needs to 
be done one time. From the bowtie_genomes folder run:

bowtie2-build select.fa select

A manifest file contains the list of files that need to be run through
the pipeline. It contains six space separated columns. Lines with a #
are ignored.

read1.fastq read2.fastq setName expectedLength sequence alignName

Within the 00_run_pipeline.py file the location to the bowtie2 index 
and the number of available threads are hardcoded.
To run the pipeline:
python 00_run_pipeline.py manifest.txt

The final mutations and counts are in the setName_mutfind.txt file.
Mutations are in the REVERSE COMPLEMENT and will to be converted to
the forward sequence. They are all relative to the provided reference
sequence and are 0-indexed. For example 0-23G: the first position is
deleted, the 23rd position is a substitution to G.
