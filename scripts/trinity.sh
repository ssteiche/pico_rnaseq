#!/usr/bin/env bash
#BSUB -J trinity
#BSUB -o logs/trinity_%J.out
#BSUB -e logs/trinity_%J.err
#BSUB -R "select[mem>150] rusage[mem=150] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source config.sh

$TRINITY \
--genome_guided_bam /mnt/alignments/temperature_expt_sorted.bam \
--max_memory 100G \
--genome_guided_max_intron 10000 \
--genome_guided_min_coverage 50 \
--genome_guided_min_reads_per_partition 50 \
--CPU 8 \
--output /mnt/trinity_assembly

