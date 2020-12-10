#!/usr/bin/env bash
#BSUB -J stringtie1[1-62]%5
#BSUB -o logs/stringtie_step1_%J.out
#BSUB -e logs/stringtie_step1_%J.err
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
	bam=/mnt/alignments/${SAMPLE}_sorted.bam
	$stringtie $bam \
	-o /mnt/stringtie/$SAMPLE.gtf \
	-p 12 \
	-G /mnt/$gtf_ref \
	-m 150 \
	-c 5
done

# need to follow up with:

# stringtie --merge
# gffcompare
# stringtie -eB

