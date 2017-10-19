#!/usr/bin/env bash
#BSUB -J run
#BSUB -o logs/run_%J.out
#BSUB -e logs/run_%J.err

bsub < code/stringtie_step1.sh
bsub -w "done(stringtie1)" < code/stringtie_merge.sh
bsub -w "done(stringtie2)" < code/stringtie_step3.sh

