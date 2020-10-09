#!/usr/bin/env bash
#BSUB -J build
#BSUB -o logs/gmap_build_%J.out
#BSUB -e logs/gsnap_build_%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q normal
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source code/config.sh

$GSNAP_DIR/gmap_build -d $SPECIES $fasta

cat $gtf_ref | $GSNAP_DIR/gtf_splicesites > $SPECIES.splicesites
cat $gtf_ref | $GSNAP_DIR/gtf_introns > $SPECIES.introns
cat $SPECIES.splicesites | $GSNAP_DIR/iit_store -o $SPECIES.splicesites.iit
cat $SPECIES.introns | $GSNAP_DIR/iit_store -o $SPECIES.introns.iit

mv $SPECIES.splicesites $MAPDIR/$SPECIES
mv $SPECIES.introns $MAPDIR/$SPECIES
mv $SPECIES.splicesites.iit $MAPDIR/$SPECIES/$SPECIES.maps
mv $SPECIES.introns.iit $MAPDIR/$SPECIES/$SPECIES.maps

