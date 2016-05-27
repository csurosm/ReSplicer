The following is a complete example of analyzing and reannotating gene structures, 
in lieu of a long fancy manual. 

The test data cover nine oomycete genomes. The analysis results are described in our paper: 

>Bocco and Cs&#369;rös "Splice sites seldom slide: intron evolution in oomycetes", Genome Biology and Evolution, 2016.

# 0. Collect the data

## 0.1 Downloads

The genome sequences and the annotations can be downloaded from the URLs below. JGI 
requires that you register with them, but others can be downloaded directly from the command line
(e.g., `wget -O albu.fna.gz ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/albugo_laibachii//dna/Albugo_laibachii.ENA1.31.dna_sm.genome.fa.gz`).
Rename the files as shown (four-letter organism code + extension). You can leave the files compressed 
(here, JGI and Broad files are uncompressed, those from Ensembl Protists stay gzipped). 

Note: P. parasitica genome has a newer assembly. We used the originally published files 
(phytophthora_parasitica_inra-310_2_supercontigs.fasta, phytophthora_parasitica_inra-310_2_transcripts.gtf),
but Broad redesigned their data download pages since. 

### albu
* albu.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/albugo_laibachii/dna/Albugo_laibachii.ENA1.31.dna_sm.genome.fa.gz
* albu.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/albugo_laibachii/Albugo_laibachii.ENA1.31.gff3.gz

### hyal
* hyal.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/hyaloperonospora_arabidopsidis/dna/Hyaloperonospora_arabidopsidis.HyaAraEmoy2_2.0.31.dna_sm.genome.fa.gz
* hyal.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/hyaloperonospora_arabidopsidis/Hyaloperonospora_arabidopsidis.HyaAraEmoy2_2.0.31.gff3.gz

### phca
* phca.fna	 http://genome.jgi.doe.gov/Phyca11/download/Phyca11_masked_genomic_scaffolds.fasta.gz
* phca.gff	 http://genome.jgi.doe.gov/Phyca11/download/Phyca11_filtered_genes.gff.gz

### phci
* phci.fna	 http://genome.jgi.doe.gov/Phyci1/download/Phyci1_AssemblyScaffolds.fasta.gz
* phci.gff	 http://genome.jgi.doe.gov/Phyci1/download/Phyci1_GeneCatalog_genes_20120612.gff.gz

### phin
* phin.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_infestans/dna/Phytophthora_infestans.ASM14294v1.31.dna_sm.genome.fa.gz
* phin.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_infestans/Phytophthora_infestans.ASM14294v1.31.gff3.gz

### phpa
* phpa.fna	 https://olive.broadinstitute.org/collections/phytophthora_parasitica_inra_310.2/downloads/scaffolds.fasta
* phpa.gtf	 https://olive.broadinstitute.org/collections/phytophthora_parasitica_inra_310.2/downloads/genes.gtf	

### phra
* phra.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_ramorum/dna/Phytophthora_ramorum.ASM14973v1.31.dna_sm.genome.fa.gz
* phra.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_ramorum/Phytophthora_ramorum.ASM14973v1.31.gff3.gz

### phso
* phso.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_sojae/dna/Phytophthora_sojae.ASM14975v1.31.dna_sm.genome.fa.gz
* phso.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_sojae/Phytophthora_sojae.ASM14975v1.31.gff3.gz

### pyul
* pyul.fna.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/pythium_ultimum/dna/Pythium_ultimum.pug.31.dna_sm.genome.fa.gz
* pyul.gff.gz	 ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/pythium_ultimum/Pythium_ultimum.pug.31.gff3.gz
	
## 0.2 Orthologs

The ortholog sets are compiled in the tarball oo_orthologs.tar.gz. Download and unpack (tar zxf oo_orthologs.tar.gz): it puts the data files in a directory called oo_orthologs.

Now your directory should look like this: 

	% ls -lph 
	total 861472
	-rw-r--r--@    1 csuros  staff   6.7K 27 May 13:19 README.md
	-rw-r--r--     1 csuros  staff    10M 27 May 12:06 albu.fna.gz
	-rw-r--r--     1 csuros  staff   1.8M 27 May 12:14 albu.gff.gz
	-rw-r--r--     1 csuros  staff    22M 27 May 12:13 hyal.fna.gz
	-rw-r--r--     1 csuros  staff   1.7M 27 May 12:13 hyal.gff.gz
	drwxr-xr-x  1919 csuros  staff    64K 27 May 00:08 oo_orthologs/
	-rw-r--r--     1 csuros  staff   3.4M 27 May 00:09 oo_orthologs.tar.gz
	-rw-r--r--@    1 csuros  staff   229B 27 May 12:21 oomycetes.tre
	-rw-r--r--@    1 csuros  staff    62M 27 May 12:10 phca.fna
	-rw-r--r--@    1 csuros  staff    12M 27 May 12:11 phca.gff
	-rw-r--r--     1 csuros  staff    75M 27 May 11:52 phci.fna
	-rw-r--r--@    1 csuros  staff    13M 27 May 12:04 phci.gff
	-rw-r--r--     1 csuros  staff    58M 27 May 11:53 phin.fna.gz
	-rw-r--r--     1 csuros  staff   4.7M 27 May 12:08 phin.gff.gz
	-rw-r--r--     1 csuros  staff    80M 10 Feb  2014 phpa.fna
	-rw-r--r--     1 csuros  staff    16M 10 Feb  2014 phpa.gtf
	-rw-r--r--     1 csuros  staff    17M 27 May 11:47 phra.fna.gz
	-rw-r--r--     1 csuros  staff   1.9M 27 May 12:08 phra.gff.gz
	-rw-r--r--     1 csuros  staff    24M 27 May 11:49 phso.fna.gz
	-rw-r--r--     1 csuros  staff   2.7M 27 May 12:07 phso.gff.gz
	-rw-r--r--     1 csuros  staff    13M 27 May 12:13 pyul.fna.gz
	-rw-r--r--     1 csuros  staff   2.0M 27 May 12:13 pyul.gff.gz


# 1. Combine the gene structure annotations

Add all the (GFF, GTF) annotations into a single file that the programs can use. The example shows that you can process multiple annotation
files if they have the same format. 

	java -cp ../ReSplicer.jar splice.extractAnnotations -prefix '|' -o annotations.txt gff3 albu=albu.fna.gz,hyal=hyal.fna.gz,phin=phin.fna.gz,phra=phra.fna.gz,phso=phso.fna.gz,pyul=pyul.fna.gz albu=albu.gff.gz,hyal=hyal.gff.gz,phin=phin.gff.gz,phra=phra.gff.gz,phso=phso.gff.gz,pyul=pyul.gff.gz
	java -cp ../ReSplicer.jar splice.extractAnnotations -prefix '|' jgi phca=phca.fna,phci=phci.fna phca=phca.gff,phci=phci.gff  >> annotations.txt 
	java -cp ../ReSplicer.jar splice.extractAnnotations -prefix '|' gtf phpa=phpa.fna phpa=phpa.gtf >> annotations.txt 

(You will get warnings about short, 0-3nt introns.)

# 2. Build the models

The example shows that you may compress the model data

	java -cp ../ReSplicer.jar splice.collectStatistics albu=albu.fna.gz,hyal=hyal.fna.gz,phca=phca.fna,phci=phci.fna,phin=phin.fna.gz,phpa=phpa.fna,phra=phra.fna.gz,phso=phso.fna.gz,pyul=pyul.fna.gz annotations.txt oo_orthologs/*.fasta | gzip - > oo.data.gz

# 3. Ancestral reconstruction

The log file is written into respliced.log, the updated gene information into a directory called respliced/.

	java -cp ../ReSplicer.jar splice.checkSites -o respliced.log -save respliced albu=albu.fna.gz,hyal=hyal.fna.gz,phca=phca.fna,phci=phci.fna,phin=phin.fna.gz,phpa=phpa.fna,phra=phra.fna.gz,phso=phso.fna.gz,pyul=pyul.fna.gz annotations.txt oo.data.gz oomycetes.tre oo_orthologs/*.fasta
	
Your directory should look like this at the end: 

	% ls -lph
	total 920848
	-rw-r--r--@    1 csuros  staff   6.7K 27 May 13:19 README.md
	-rw-r--r--     1 csuros  staff    10M 27 May 12:06 albu.fna.gz
	-rw-r--r--     1 csuros  staff   1.8M 27 May 12:14 albu.gff.gz
	-rw-r--r--     1 csuros  staff    20M 27 May 12:52 annotations.txt
	-rw-r--r--     1 csuros  staff    22M 27 May 12:13 hyal.fna.gz
	-rw-r--r--     1 csuros  staff   1.7M 27 May 12:13 hyal.gff.gz
	-rw-r--r--     1 csuros  staff   3.1M 27 May 13:00 oo.data.gz
	drwxr-xr-x  1919 csuros  staff    64K 27 May 00:08 oo_orthologs/
	-rw-r--r--     1 csuros  staff   3.4M 27 May 00:09 oo_orthologs.tar.gz
	-rw-r--r--@    1 csuros  staff   229B 27 May 12:21 oomycetes.tre
	-rw-r--r--@    1 csuros  staff    62M 27 May 12:10 phca.fna
	-rw-r--r--@    1 csuros  staff    12M 27 May 12:11 phca.gff
	-rw-r--r--     1 csuros  staff    75M 27 May 11:52 phci.fna
	-rw-r--r--@    1 csuros  staff    13M 27 May 12:04 phci.gff
	-rw-r--r--     1 csuros  staff    58M 27 May 11:53 phin.fna.gz
	-rw-r--r--     1 csuros  staff   4.7M 27 May 12:08 phin.gff.gz
	-rw-r--r--     1 csuros  staff    80M 10 Feb  2014 phpa.fna
	-rw-r--r--     1 csuros  staff    16M 10 Feb  2014 phpa.gtf
	-rw-r--r--     1 csuros  staff    17M 27 May 11:47 phra.fna.gz
	-rw-r--r--     1 csuros  staff   1.9M 27 May 12:08 phra.gff.gz
	-rw-r--r--     1 csuros  staff    24M 27 May 11:49 phso.fna.gz
	-rw-r--r--     1 csuros  staff   2.7M 27 May 12:07 phso.gff.gz
	-rw-r--r--     1 csuros  staff    13M 27 May 12:13 pyul.fna.gz
	-rw-r--r--     1 csuros  staff   2.0M 27 May 12:13 pyul.gff.gz
	drwxr-xr-x  1322 csuros  staff    44K 27 May 13:09 respliced/
	-rw-r--r--     1 csuros  staff   6.0M 27 May 13:09 respliced.log


# 4. Questions / problems 

Send me an e-mail: csurosm@gmail.com

Miklós Cs&#369;rös
