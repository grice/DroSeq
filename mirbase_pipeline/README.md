# Droseq mirbase pipeline

REQUIREMENTS python3

Command line programs:
 - sickle
 - bowtie2
 - samtools

Python modules:
 - pysam

Bowtie2 database containing pri-miRNA sequences in reverse complement
based on the library construction strategy. This step only needs to 
be done one time. From the bowtie_genomes folder run:

bowtie2-build combo_rc_padded.fa combo_rc_padded

A manifest file contains the list of files that need to be run through
the pipeline. It contains three space separated columns. Lines with a #
are ignored.

fastq_r1.fastq.gz fastq_r2.fastq.gz sampleName

Within the 00_run_pipeline.py file the location to the bowtie2 index 
and the number of available threads are hardcoded.
To run the pipeline:
python 00_run_pipeline.py manifest.txt
