

# ReSplicer

The package implements our analysis methods for exon-intron structure evolution. The program evaluates splice site homology using multiple, annotated genomic sequences and aligned orthologous protein-coding genes. Ancestral gene structure reconstruction considers intron gain and loss, as well as intron boundary sliding. The program may be used also to recalculate exon boundaries based on alignment evidence.    
The software is free under the terms of the MIT License (see `LICENSE.md`).

In order to use it, you need: 

* a set of aligned orthologs  
* a phylogeny
* genome sequences

Download ReSplicer.jar. The test directory contains a complete example for 9 oomycete genomes. 


0. Collect data
===============
For every organism: get the genomes (Fasta with chromosome sequences), the gene annotations (GFF, GenePred, GTF formats work),and the proteins. Select your ortholog groups (e.g., use [OrthoMCL](http://www.orthomcl.org/)): they must contain exactly one gene from each genome, and align the protein sequences (e.g., use [Muscle](http://www.drive5.com/muscle/)). Compute a phylogeny for the genomes (branch lengths are immaterial): use a simple codename for every taxon. 

Below, `JAVA` means your java engine with appropriate options (like `java -Xmx2048M` for 2G memory). All programs can be run without arguments for more information about their arguments (e.g., `JAVA -cp ReSplicer.jar splice.extractAnnotations`). 

1. Extract gene structure annotations
=====================================
Collect all gene structure annotations into a single file.

`JAVA -cp ReSplicer.jar splice.extractAnnotations [options] format org=genome.fasta  org=genes.ext`

* format may be `gff3`, `gtf`, `gp`, or `jgi` 
* org is the code for your organism
* genes.ext is the annotation file in the given format 

Example: `java -cp ReSplicer.jar splice.extractAnnotations gff3 albu=genome-albu.fa albu=Albugo_laibachii/Albugo_laibachii.ENA1.26.gff3 >> annotations.txt`

2. Build the models 
===================
Using the alignments and the gene structure annotations, the program calculates statistics about codon substitutions and gene structures.

`JAVA -cp ReSplicer.jar splice.collectStatistics [options] org1=g1.fa,org2=g2.fa,... annotations.txt ali1.fa ali2.fa ... `

* org=g.fa specifies the Fasta for the genome of org
* annotations.txt is the collected gene structure info from Step 1
* ali1, ali2, ... are your aligned orthologs in Fasta format

Example: `java -Xmx2048M-cp ReSplicer.jar splice.collectStatistics -o model.data phra=genome-phra.fa.gz,phso=genome-phso.fa.gz,pyul=genome-pyul.fa.gz,albu=genome-albu.fa,hyal=genome-hyal.fa,phca=genome-phca.fa,phin=genome-phin.fa.gz,phpa=genome-phpa.fa.gz,phci=genome-phci.fa.gz annotations.txt orthologs/ali*.faa`

3. Do the reconstruction
========================
`JAVA -cp ReSplicer.jar splice.checkSites [options] org1=g1.fa,org2=g2.fa,... annotations.txt tree model.data ali1.fa ali2.fa ... `

* tree is a [Newick-format](http://evolution.genetics.washington.edu/phylip/newicktree.html) tree
* options include `-save` to write reannotated genes into a separate directory

Example: `java -Xmx2048M -cp ReSplicer.jar splice.checkSites -loss 5 -gain 12 -annot 1 -shift 5 -save savedir  phra=genome-phra.fa.gz,phso=genome-phso.fa.gz,pyul=genome-pyul.fa.gz,albu=genome-albu.fa,hyal=genome-hyal.fa,phca=genome-phca.fa,phin=genome-phin.fa.gz,phpa=genome-phpa.fa.gz,phci=genome-phci.fa.gz annotations.txt model.data oomycetes.tre orthologs/ali*.faa`

