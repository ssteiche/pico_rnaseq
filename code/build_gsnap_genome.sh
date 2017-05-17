#!/usr/bin/env bash
#BSUB -J build
#BSUB -o logs/gmap_build_%J.out
#BSUB -e logs/gsnap_build_%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q normal
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

genome_name=se107_orig
GSNAP_DIR=$HOME/bin/gsnap_2017-05-08/bin
MAPDIR=$HOME/gmapdb
fasta=$HOME/Projects/llaurens/se107/assembly/SE107.scaffolds.fa
#gtf_ref=$HOME/Projects/llaurens/se107/annotation/SE107_augustus_gene_predictions.gtf
gtf_ref=genome/se107_augustus_exons.gtf

$GSNAP_DIR/gmap_build -d $genome_name $fasta

cat $gtf_ref | $GSNAP_DIR/gtf_splicesites > $genome_name.splicesites
cat $gtf_ref | $GSNAP_DIR/gtf_introns > $genome_name.introns
cat $genome_name.splicesites | $GSNAP_DIR/iit_store -o $genome_name.splicesites.iit
cat $genome_name.introns | $GSNAP_DIR/iit_store -o $genome_name.introns.iit

mv $genome_name.splicesites $MAPDIR/$genome_name
mv $genome_name.introns $MAPDIR/$genome_name
mv $genome_name.splicesites.iit $MAPDIR/$genome_name/$genome_name.maps
mv $genome_name.introns.iit $MAPDIR/$genome_name/$genome_name.maps

