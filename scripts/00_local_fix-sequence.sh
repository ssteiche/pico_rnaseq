#!/usr/bin/env bash
#BSUB -J gsnap[1-62]
#BSUB -o logs/preprocess_%J.out
#BSUB -e logs/preprocess_%J.err
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

#62 samples

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

# old script pointed to code/config.sh
source config.sh

for SAMPLE in ${SAMPLES[@]}
do
	gunzip -c ${fastq_folder}/FilteredReads/${SAMPLE}_sequence.txt.gz > ${fastq_folder}/${SAMPLE}.fastq
	gzip ${fastq_folder}/${SAMPLE}.fastq
done
