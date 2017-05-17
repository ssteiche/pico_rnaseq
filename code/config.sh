#!/usr/bin/env bash

alignment_dir=alignments
GSNAP_DIR=$HOME/bin/gsnap_2017-05-08/bin
MAPDIR=$HOME/gmapdb
MAPPER=$GSNAP_DIR/gsnap
fastq_folder=fastq
SPECIES=se107_orig
SNPS=""
MAX_MM="--max-mismatches=0.06"
PLATFORM=illumina
KNOWN_SPLICE="--use-splicing=se107_orig.splicesites"

if [[ ! -d $alignment_dir ]]; then
        mkdir -p $alignment_dir
fi
if [[ ! -d logs ]]; then
        mkdir logs
fi

SAMPLES=(
0hr_RT_P1
0hr_RT_P2
1hr_C_P1
1hr_C_P2
1hr_H_P1
1hr_H_P2
1hr_RT_P1
1hr_RT_P2
1hr_S1_Con
1hr_S1_FD
1hr_S1_Rot
1hr_S2_Con
1hr_S2_FD
1hr_S2_Rot
24hr_C_P1
24hr_C_P2
24hr_H_P1
24hr_H_P2
24hr_RT_P1
24hr_RT_P2
24hr_S1_Con
24hr_S1_FD
24hr_S1_Rot
24hr_S2_Con
24hr_S2_FD
24hr_S2_Rot
3hr_C_P1
3hr_C_P2
3hr_H_P1
3hr_H_P2
3hr_RT_P1
3hr_RT_P2
3hr_S1_Con
3hr_S1_FD
3hr_S1_Rot
3hr_S2_Con
3hr_S2_FD
3hr_S2_Rot
72hr_C_P1
72hr_C_P2
72hr_H_P1
72hr_H_P2
72hr_RT_P1
72hr_RT_P2
72hr_S1_Con
72hr_S1_FD
72hr_S1_Rot
72hr_S2_Con
72hr_S2_FD
72hr_S2_Rot
9hr_C_P1
9hr_C_P2
9hr_H_P1
9hr_H_P2
9hr_RT_P1
9hr_RT_P2
9hr_S1_Con
9hr_S1_FD
9hr_S1_Rot
9hr_S2_Con
9hr_S2_FD
9hr_S2_Rot
)



