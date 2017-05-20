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

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=fastq/${SAMPLE}.fastq.gz
bam=alignments/${sample}_sorted.bam


stringtie $bam \
-o stringtie/$sample \
-p 12 \
-G ??path-to-gtf?? \




stringtie --merge

gffcompare


stringtie -eB






