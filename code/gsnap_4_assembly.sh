#!/usr/bin/env bash
#BSUB -J gsnap
#BSUB -o logs/gsnap_%J.out
#BSUB -e logs/gsnap_%J.err
#BSUB -R "select[mem>48] rusage[mem=48] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source code/config.sh

SAMPLE=temperature_expt
R1=fastq/*_P[123].fastq.gz
R2=""

# like previous run, but force --single end so that multiple fastqs are treated as single reads.
# also keep failed alignments so trinity can use them

#$MAPPER -t 12 -d $SPECIES --gunzip --force-single-end $KNOWN_SPLICE $MAX_MM --speed 2 --use-shared-memory 0 -n 1 --novelsplicing 1 -i 2 --read-group-id=$SAMPLE --read-group-name=$SAMPLE --read-group-library=$SAMPLE --read-group-platform=illumina -A sam $R1 $R2 > alignments/$SAMPLE.sam

#samtools view -hbS alignments/$SAMPLE.sam > alignments/$SAMPLE.bam

samtools sort -@ 12 -m 4G alignments/$SAMPLE.bam alignments/${SAMPLE}_sorted

