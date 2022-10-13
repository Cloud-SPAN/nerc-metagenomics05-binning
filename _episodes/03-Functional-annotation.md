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
cd ~/cs_course/analysis/
mkdir prokka
cd prokka
prokka -h
~~~
{: .bash}

From this we can build our basic command.
~~~
prokka --outdir mydir --prefix mygenome contigs.fa
~~~
{: .output}

Prokka produces multiple different file types which you can see in the table below. We are mainly interested in `.faa` and `.txt` but many of the other files are useful for submission to different databases.

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

We are going to use Prokka to initially annotate one MAG.

In the previous episode we saw we'd produced 6 MAGs of differing quality.
In our example, we are going to start with the MAG `bin.6.fa` as it was assessed by CheckM as being 99.45% complete and 0% contaminated.

~~~
prokka --outdir bin.6 --prefix bin.6 ../binning/pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/bin.6.fa
~~~
{: .bash}

This should take around a minute on the instance so we will not be running the command in the background.

When you initially run the command you should see similar to the following.
~~~
[15:45:37] This is prokka 1.14.6
[15:45:37] Written by Torsten Seemann <torsten.seemann@gmail.com>
[15:45:37] Homepage is https://github.com/tseemann/prokka
[15:45:37] Local time is Thu Oct 13 15:45:37 2022
[15:45:37] You are csuser
[15:45:37] Operating system is linux
[15:45:37] You have BioPerl 1.7.8
Argument "1.7.8" isn't numeric in numeric lt (<) at /home/csuser/.miniconda3/envs/prokka/bin/prokka line 259.
[15:45:37] System has 8 cores.
[15:45:37] Will use maximum of 8 cores.
[15:45:37] Annotating as >>> Bacteria <<<
[15:45:37] Generating locus_tag from 'bin.6.fa' contents.
~~~
{: .output}
And you should see the following when the command has finished:

~~~
[15:45:59] Output files:
[15:45:59] bin.6/bin.6.txt
[15:45:59] bin.6/bin.6.log
[15:45:59] bin.6/bin.6.tsv
[15:45:59] bin.6/bin.6.fsa
[15:45:59] bin.6/bin.6.fna
[15:45:59] bin.6/bin.6.sqn
[15:45:59] bin.6/bin.6.faa
[15:45:59] bin.6/bin.6.gbk
[15:45:59] bin.6/bin.6.ffn
[15:45:59] bin.6/bin.6.err
[15:45:59] bin.6/bin.6.tbl
[15:45:59] bin.6/bin.6.gff
[15:45:59] Annotation finished successfully.
[15:45:59] Walltime used: 0.37 minutes
[15:45:59] If you use this result please cite the Prokka paper:
[15:45:59] Seemann T (2014) Prokka: rapid prokaryotic genome annotation. Bioinformatics. 30(14):2068-9.
[15:45:59] Type 'prokka --citation' for more details.
[15:45:59] Thank you, come again.
~~~
{: .output}

If we navigate into the output file created by Prokka we can then list and see that Prokka has generated

~~~
cd bin.6
ls
~~~
{: .bash}

~~~
bin.6.err  bin.6.faa  bin.6.ffn  bin.6.fna  bin.6.fsa  bin.6.gbk  bin.6.gff  bin.6.log  bin.6.sqn  bin.6.tbl  bin.6.tsv  bin.6.txt
~~~
{: .output}

If we refer back to the Prokka documentation, we can remind ourselves what each of these files are for.

| Extension |                                                                                              Description                                                                                              |                                                                                                                             What does this mean?                                                                                                                            |
|:---------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| .gff      | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV.                                                                 | This is a standard bioinformatics format used to collate genome features of an annotation. As it suggests you can view this in IGV. You will need this file later on for CIRCOS.                                                                                            |
| .gbk      | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence.                         | A file that is created in order to upload your sequence to Genbank. We don’t want to do this at this point so you can just ignore this.                                                                                                                                     |
| .fna      | Nucleotide FASTA file of the input contig sequences.                                                                                                                                                  | FASTA file of the sequences you passed to prokka (this should be basically identical to the MAG FASTA file)                                                                                                                                                                 |
| .faa      | Protein FASTA file of the translated CDS sequences.                                                                                                                                                   | The protein sequences (i.e. amino acid code) of all identified CDS (i.e. open reading frames/ genes)                                                                                                                                                                        |
| .ffn      | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)                                                                                                            | Similar to the *.faa file but it is nucleotide (i.e. ACTG) rather than protein and this one contains the rRNA (16s etc), tRNA etc sequences identified.                                                                                                                     |
| .sqn      | An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc.                                                          | Only relevant if this was going to be directly uploaded to Genbank, which it isn’t.                                                                                                                                                                                         |
| .fsa      | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. | An intermediate file needed for one of the Prokka steps.                                                                                                                                                                                                                    |
| .tbl      | Feature Table file, used by "tbl2asn" to create the .sqn file.                                                                                                                                        | Another intermediate file needed for one of the Prokka steps.                                                                                                                                                                                                               |
| .err      | Unacceptable annotations - the NCBI discrepancy report.                                                                                                                                               | This is only relevant if you wanted to upload the prokka output to NCBI (which we don’t…)                                                                                                                                                                                   |
| .log      | Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled.                                                      | This file outlines all the steps Prokka has gone through. Useful to know what Prokka has actually done and also identified.                                                                                                                                                 |
| .txt      | Statistics relating to the annotated features found.                                                                                                                                                  | Basic output stats (i.e. counts of what have been generated) for the genome.                                                                                                                                                                                                |
| .tsv      | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product                                                                                                                 | This is a tsv file (think a simple version of an excel spreadsheet where different cells are indicated by tabs and new lines) containing all the annotations prokka has generated. This might be useful downstream if, for example, you’re looking for a specific protein.  |

Do we want to identify the 16S and then use this to make a tree?


How we can use kofamscan to generate KEGG ids

How to use KEGG decoder to plot the pathways in your assembly




{% include links.md %}
