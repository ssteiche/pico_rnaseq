#!/usr/bin/env bash
#BSUB -J build
#BSUB -o logs/gmap_build_%J.out
#BSUB -e logs/gsnap_build_%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q normal
#BSUB -P llaurens

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source config.sh

singularity exec --bind $PWD:$PWD $HOME/containers/prna.sif gmap_build -d se107_mik mikado/SE107.scaffolds.fa -D gmapdb

# first convert to gtf, because util cannot parse gff3
/usr/local/gffread/gffread mikado.loci.gff3 -T -o mikado.loci.gtf
cat mikado.loci.gtf | singularity exec --bind $PWD:$PWD $HOME/containers/prna.sif gtf_splicesites > se107_mik.splicesites
cat mikado.loci.gtf | singularity exec --bind $PWD:$PWD $HOME/containers/prna.sif gtf_introns > se107_mik.introns

cat se107_mik.splicesites | singularity exec --bind $PWD:$PWD $HOME/containers/prna.sif iit_store -o se107_mik.splicesites.iit
cat se107_mik.introns | singularity exec --bind $PWD:$PWD $HOME/containers/prna.sif iit_store -o se107_mik.introns.iit

mv se107_mik.splicesites gmapdb/se107_mik
mv se107_mik.introns gmapdb/se107_mik
mv se107_mik.splicesites.iit gmapdb/se107_mik/se107_mik.maps
mv se107_mik.introns.iit gmapdb/se107_mik/se107_mik.maps

