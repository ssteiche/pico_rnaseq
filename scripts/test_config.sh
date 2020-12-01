#!/usr/bin/env bash
#BSUB -J gsnap[1-62]
#BSUB -o logs/preprocess_%J.out
#BSUB -e logs/preprocess_%J.err
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

#set -o nounset -o pipefail -o errexit -x

source config.sh

echo $fastq_folder

# Using parameter expansion on an array input to iterate commands. Recall that [@] expands the whole array
for SAMPLE in ${SAMPLES[@]}
do
	printf "${SAMPLE}_seq.fasta\n"
done

#print $SAMPLES
