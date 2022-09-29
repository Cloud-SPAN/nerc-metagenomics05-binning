---
title: "Functional annotation"
teaching: 40
exercises: 10
questions:
- "How can we add functional annotation to our bins?"
- "How can we identify what pathways these are involved with?"
objectives:
- "What is functional annotation"
- "How to use prokka for functional annotation"

keypoints:
- "Functional annotation allows us to look at the metabolic capacity of a metagenome"  
- "Prokka can be used to predict genes in our assembly"
- "Kofam scan can be used to identify KEGG IDs and enzyme numbers"
math: true
---

## What is functional annotation?


## How we perform functional annotation?

We will be annotating each of our MAGs using [Prokka](https://github.com/tseemann/prokka) for rapid prokaryotic genome annotation on the command line.

> ## Software choices
> We are using prokka here as it is still the software most commonly used. However, the program is no longer being updated. One recent alternative that is being actively developed is [Bakta](https://github.com/oschwengers/bakta).
{: .callout}

Prokka identifies candidate genes in a iterative process. First using Prodigal (another command line tool that prokka uses in the pipeline) to find candiate genes and then they are then compared against databases of known protein sequences in order to determine their function. You can read more about Prokka in the corresponding paper [Seeman, 2014](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).

Prokka has been pre-installed on our instance and we can access the help documentation using
~~~
prokka -h
~~~
{: .bash}

From this we can see build our command.
~~~
prokka --outdir mydir --prefix mygenome contigs.fa
~~~
{: .bash}

We need to run prokka for each of our six MAGs. 

| Suffix | Description of file contents                       |
|--------|----------------------------------------------------|
| .fna   | FASTA file of original input contigs (nucleotide)  |
| .faa   | FASTA file of translated coding genes (protein)    |
| .ffn   | FASTA file of all genomic features (nucleotide)    |
| .fsa   | Contig sequences for submission (nucleotide)       |
| .tbl   | Feature table for submission                       |
| .sqn   | Sequin editable file for submission                |
| .gbk   | Genbank file containing sequences and annotations  |
| .gff   | GFF v3 file containing sequences and annotations   |
| .log   | Log file of Prokka processing output               |
| .txt   | Annotation summary statistics                      |


Do we want to identify the 16S and then use this to make a tree?


How we can use kofamscan to generate KEGG ids

How to use KEGG decoder to plot the pathways in your assembly




{% include links.md %}
