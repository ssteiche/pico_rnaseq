#!/usr/bin/env bash
#BSUB -J gsnap
#BSUB -o logs/gsnap_%J.out
#BSUB -e logs/gsnap_%J.err
#BSUB -R "select[mem>48] rusage[mem=48] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source config.sh

SAMPLE=temperature_expt
#R1=/mnt/*_P[123].fastq.gz
R2=""

# like previous run, but force --single end so that multiple fastqs are treated as single reads.
# also keep failed alignments so trinity can use them

$MAPPER -D $HOME/projects/gmapdb -d $SPECIES -t 8 --gunzip --force-single-end $KNOWN_SPLICE $MAX_MM --use-shared-memory 0 -n 1 --novelsplicing 1 -i 2 --read-group-id=$SAMPLE --read-group-name=$SAMPLE --read-group-library=$SAMPLE --read-group-platform=illumina -A sam $TEMP_SAMPLES $R2 > $alignment_dir/$SAMPLE.sam

$SAMTOOLS view -hbS /mnt/alignments/$SAMPLE.sam > $alignment_dir/$SAMPLE.bam

$SAMTOOLS sort -@ 8 -m 4G /mnt/alignments/$SAMPLE.bam -o /mnt/alignments/${SAMPLE}_sorted

rm $alignment_dir/$SAMPLE.sam
rm $alignment_dir/$SAMPLE.bam

