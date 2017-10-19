
# Notes on using Trinity

Installed at `$HOME/bin/Trinity-v2.4.0`

1. Run GSNAP on using all temp experiments

2. Sort using samtools

3. Run Trinity

```{bash}
 Trinity \
  --genome_guided_bam rnaseq.coordSorted.bam \
  --genome_guided_max_intron 10000 \
  --max_memory 40G \
  --CPU 12 
```




