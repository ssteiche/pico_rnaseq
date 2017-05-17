#!/usr/bin/env bash
#BSUB -J gsnap[1-62]%5
#BSUB -o logs/gsnap_%J.out
#BSUB -e logs/gsnap_%J.err
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

#62 samples

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source code/config.sh

SAMPLE=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
#R1=fastq/FilteredReads/${SAMPLE}_sequence.txt.gz
R1=fastq/${SAMPLE}.fastq.gz
R2=""

$MAPPER -m 5 --gunzip -d $SPECIES $KNOWN_SPLICE $MAX_MM -i2 -n1 -Q --nofails -B5 -t12 --read-group-id=$SAMPLE --read-group-name=$SAMPLE --read-group-library=$SAMPLE --read-group-platform=illumina -A sam $R1 $R2 > alignments/$SAMPLE.sam

samtools view -hbS alignments/$SAMPLE.sam > alignments/$SAMPLE.bam

