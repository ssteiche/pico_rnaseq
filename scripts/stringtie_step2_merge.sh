#!/usr/bin/env bash
#BSUB -J stringtie2
#BSUB -o logs/stringtie_merge_%J.out
#BSUB -e logs/stringtie_merge_%J.err
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 12
#BSUB -P llaurens

#62 samples

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source config.sh

ls $stringtie_dir/*.gtf > $stringtie_dir/merged_list.txt

$stringtie --merge $stringtie_dir/merged_list.txt \
-o /mnt/f/genomics/scen_rnaseq/genome/stringtie_merge.gtf \
-p 12 \
-G /mnt/f/genomics/scen_rnaseq/$gtf_ref \
-m 150 \
-c 5 \
-F 5

### Not currently installed in container. Add then re-do
#$GFFCOMPARE -r /mnt/f/genomics/scen_rnaseq/$gtf_ref -G -o compare /mnt/f/genomics/scen_rnaseq/genome/stringtie_merge.gtf


# followup with
# stringtie -eB

