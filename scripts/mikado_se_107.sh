# Rough script for testing mikado transcript integration on se_107 dataset
# Minimally need two gtf-like files to compare
# Also request a set of high quality junctions, produced below by Portcullis

cd /mnt/f/genomics/scen_rnaseq/portcullis

# Inputs are the original scaffold assembly and the combined alignments of the temperature experiment RNA-seq
portcullis full -t 8 ../genome/SE107.scaffolds.fa ../alignments/temperature_expt_sorted.bam

cd ../mikado
mikado configure --list list.txt --reference SE107.scaffolds.fa --mode permissive --scoring plant.yaml --copy-scoring plant.yaml --junctions ../portcullis/portcullis_out/2-junc/portcullis_all.junctions.bed test1_configuration.yaml

mikado prepare --json-conf test1_configuration.yaml

# BLAST against UNIPROT plant references. Obtain by visiting uniprot an searching for "taxonomy:viridiplantae AND reviewed:yes", then download all results
gzip -dc uniprot_sp_plants_122120.fasta.gz > uniprot_sprot_plants.fasta
makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log
blastx -max_target_seqs 5 -num_threads 8 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz
# Time consumed on desktop
#real    101m48.998s
#user    584m0.549s
#sys     12m29.173s

# Use transdecoder to identify ORFs in SE genome assembly/annotation
# had to install a perl module to run
cpan URI::Escape
# Use the gtf file-based instructions here: https://github.com/TransDecoder/TransDecoder/wiki
#/usr/local/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl ../genome/se107_augustus_minimal.gtf ../genome/SE107.scaffolds.fa > td_se_transcripts.fasta
#/usr/local/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl ../genome/se107_augustus_minimal.gtf > td_se_transcripts.gff3
#TransDecoder.LongOrfs -t td_se_transcripts.fasta
# Skipped the homology prediction step
#TransDecoder.Predict -t td_se_transcripts.fasta

########### The above ORF prediction was done on the wrong file. mikado serialise requires ORFs constructed from the mikado_prepared.fasta
TransDecoder.LongOrfs -t mikado_prepared.fasta
# Skipped the homology prediction step
TransDecoder.Predict -t mikado_prepared.fasta

# The above prodiced a number of files, including the .bed used for orfs argument below
mikado serialise --json-conf test1_configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed --blast_targets uniprot_sprot_plants.fasta

# pick loci
mikado pick --json-conf test1_configuration.yaml --subloci-out mikado.subloci.gff3

# Now compare the resulting transcripts to the original
# First index the reference transcript file for quicker eval
mikado compare -r se107_augustus_minimal.gtf --index

mikado compare -r se107_augustus_minimal.gtf -p mikado_prepared.gtf -o compare_input -l compare_input.log
mikado compare -r se107_augustus_minimal.gtf -p mikado.subloci.gff3 -o compare_subloci -l compare_subloci.log
mikado compare -r se107_augustus_minimal.gtf -p mikado.loci.gff3 -o compare -l compare.log

##################################################################################################
# Follow up with alignment of reads providing snap with new mikado transcript annotations
##################################################################################################

time ./02_mik_loc_gsnap.sh
#real    681m12.104s
#user    1980m8.791s
#sys     355m31.342s

# Now that sequences are aligned again we will calculate read coverage tables with StinrTie -eB
# Note that this skips the generation of stringtie gtfs and merging. Could consider doing, but seems redundant to mikado

time ./mik_stringtie_step3.sh
