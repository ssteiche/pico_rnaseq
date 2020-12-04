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

source config.sh

for SAMPLE in ${SAMPLES[@]}
#R1=fastq/FilteredReads/${SAMPLE}_sequence.txt.gz
do
	R1=/mnt/${SAMPLE}.fastq.gz
	R2=""

	$MAPPER -D $HOME/projects/gmapdb -d $SPECIES -t 8 --gunzip $KNOWN_SPLICE $MAX_MM --use-shared-memory 0 -n 1 --novelsplicing 1 -i 2 --nofails --read-group-id=$SPECIES --read-group-name=$SPECIES --read-group-library=$SPECIES --read-group-platform=illumina -A sam $R1 $R2 > $alignment_dir/$SAMPLE.sam
	
	$SAMTOOLS view -hbS /mnt/alignments/$SAMPLE.sam > $alignment_dir/$SAMPLE.bam

	$SAMTOOLS sort /mnt/alignments/$SAMPLE.bam -@ 8 -o /mnt/alignments/${SAMPLE}_sorted.bam

	rm $alignment_dir/$SAMPLE.sam
	rm $alignment_dir/$SAMPLE.bam
done
