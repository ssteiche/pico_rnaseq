#!/usr/bin/env bash

source config.sh

# Generate CDS transcripts fasta using mikado coordinates
gffread -w $fastq_folder/genome/mik_transcripts.fa \
        -g $fastq_folder/genome/SE107.scaffolds.fa $fastq_folder/mikado.loci.gff3

salmon index -t $fastq_folder/genome/mik_transcripts.fa -i $fastq_folder/genome/mik_se107_index
