#!/usr/bin/env bash

source config.sh

for SAMPLE in ${SAMPLES[@]}
do
	echo "processing sample ${SAMPLE}"
	salmon quant -i $fastq_folder/genome/mik_se107_index -l A \
		-r $fastq_folder/${SAMPLE}.fastq.gz \
		-p 8 --validateMappings -o $fastq_folder/salmon_quants/${SAMPLE}_quant
done
