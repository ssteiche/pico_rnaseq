
# Code for processing RNAseq data for SE107 *Desmodesmus armatus*


## Quick Start

0. Setup the `config.sh` file with your settings and fix sequence reads
1. Build a genome index for SE107 for your favorite aligner
2. Align reads to the SE107 genome
3. Assemble transcripts with StringTie allowing for novel isoforms
4. Merge isoforms for each sample into a master transcriptome reference
5. Re-assemble transcripts using the new master transcriptome reference
6. Quantify transcripts with Ballgown








## Dependencies

- GSNAP
- StringTie
- samtools
- StringTie
- Ballgown





### config.sh

The config file is where all the project specific parameters and sample names
should go. The other scripts should be as abstract as possible for reuse. 


## 01\_build-gsnap-genome.sh

This script builds the necessary genome index for GSNAP.


## 02\_gsnap.sh

This script aligns reads to the SE107 genome with GSNAP


## 02\_stringtie.sh

This script is a wrapper for the following three scripts. It queues the scripts in LSF so that they are run sequentially for each sample 


### stringtie\_step1.sh

This script assembles transcripts allowing for novel isoform assembly

### stringtie\_step2\_merge.sh

This script merges the known and novel isoforms for all samples into a master transcriptome reference

### stringtie\_step3.sh

This script runs StringTie again but using the master transcriptome reference.


## gsnap\_4\_assembly.sh

This script aligns reads to the SE107 using GSNAP, but in a manner that's friendly to Trinity

## trinity.sh

This script runs Trinity which does a de-novo genome assembly.

