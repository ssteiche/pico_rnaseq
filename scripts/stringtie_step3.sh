#!/usr/bin/env bash
#BSUB -J stringtie3[1-62]%5
#BSUB -o logs/stringtie_step3_%J.out
#BSUB -e logs/stringtie_step3_%J.err
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

#62 samples

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source config.sh

#sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
#bam=alignments/${sample}_sorted.bam

for SAMPLE in ${SAMPLES[@]}
do
	bam=/mnt/f/genomics/scen_rnaseq/alignments/${SAMPLE}_sorted.bam
	$stringtie -eB \
	-p 12 \
	-G /mnt/f/genomics/scen_rnaseq/genome/stringtie_merge.gtf \
	-m 150 \
	-c 5 \
	-o /mnt/f/genomics/scen_rnaseq/ballgown/$SAMPLE/${SAMPLE}_step3.gtf \
	$bam
done

