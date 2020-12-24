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

	$MAPPER -D $mik_MAPDIR -d $mik_SPECIES -t 8 --gunzip $mik_KNOWN_SPLICE $MAX_MM --use-shared-memory 0 -n 1 --novelsplicing 1 -i 2 --nofails --read-group-id=$mik_SPECIES --read-group-name=$mik_SPECIES --read-group-library=$mik_SPECIES --read-group-platform=illumina -A sam $R1 $R2 > $mik_alignment_dir/$SAMPLE.sam
	
	$SAMTOOLS view -hbS /mnt/mikado_alignments/$SAMPLE.sam > $mik_alignment_dir/$SAMPLE.bam

	$SAMTOOLS sort /mnt/mikado_alignments/$SAMPLE.bam -@ 8 -o /mnt/mikado_alignments/${SAMPLE}_sorted.bam

	rm $mik_alignment_dir/$SAMPLE.sam
	rm $mik_alignment_dir/$SAMPLE.bam
done
