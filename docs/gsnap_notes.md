
# GNSAP Notes

Here are the args to pass when compiling GSNAP

```{bash}
./configure --prefix=$HOME/bin/gsnap_2017-05-08 --with-gmapdb=$HOME/gmapdb MAX_STACK_READLENGTH=250
make install --prefix=$HOME/bin/gsnap_2017-05-08
```

From the README on building the genome

```{bash}
gmap_build -d <genome> [-k <kmer size>] -C <fasta_files...>
```

According to the GSNAP manual it is possible to build an genome reference for scaffolds if the positions are known

    5a.  Contigs are mapped to chromosomes
    ======================================
    
    If your FASTA entries consist of contigs, each of which has a mapping
    to a chromosomal region in the header, you may add the -C (or
    --contigs-are-mapped) flag to gmap_build, like this
    
        gmap_build -d <genome> -C <fasta_file>...
    
    Then gmap_build will try to parse a chromosomal region from each
    header.  The program knows how to parse the following patterns:
    
        chr=1:217281..257582  [may insert spaces around '=', or omit '=' character]
        chr=1                 [may insert spaces around '=', or omit '=' character]
        chromosome 1                                                         [NCBI .mfa format]
        chromosome:NCBI35:22:1:49554710:1                                    [Ensembl format]
        /chromosome=2                                                        [Celera format]
        /chromosome=2 /alignment=(88840247-88864134) /orientation=rev        [Celera format]
        chr1:217281..257582
        chr1                  [may insert spaces after 'chr']
    
    If only the chromosome is specified, without coordinates, the program
    will assign its own chromosomal coordinates by concatenating the
    contigs within each chromosome.  If gmap_build cannot figure out the
    chromosome, it will assign it to chromosome "NA".

