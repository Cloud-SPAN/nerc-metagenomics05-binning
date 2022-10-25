---
title: "Metagenome Binning"
teaching: 50
exercises: 50
questions:
- "How can we obtain the original genomes from a metagenome?"
objectives:
- "Obtain Metagenome-Assembled Genomes (MAGs) from the metagenomic assembly."
- "Understand that there are multiple methods that can be used to perform binning"  
keypoints:
- "Metagenome-Assembled Genomes (MAGs) sometimes are obtained from curated contigs grouped into bins."
- "Use MetaBAT2 to assign the contigs to bins of different taxa."
- "Other programmes are available that are generating other bins, and these can be rationalised using tools such as DAStools"
---

## Metagenomic binning
Now we have an assembly we can start to separate out the individual genomes using a process called binning. This will allow us to analyse each of the species inside our sample individually. We call these genomes metagenome-assembled genomes (MAGs). 

The assembled contigs that make up the metagenome will be assigned to different "bins" (FASTA files that contain certain contigs). Ideally, each bin will correspond to only one original genome (a MAG). 

As we covered in the [assembly section](https://cloud-span.github.io/metagenomics01-qc-assembly/03-assembly/index.html), assembling the pieces of a metagenome is more difficult compared to a single genome assembly. Most assemblers are not able to reconstruct complete genomes for the organisms that are represented in the metagenome. As a result each organism wil be represented by multiple contigs following assembly and polishing. This means that we need to be able to separate these contigs so we can identify which belong to each organism in our metagenome. This is where binning comes in.

<a href="{{ page.root }}/fig/03-05-01.png">
  <img src="{{ page.root }}/fig/03-05-01.png" width="435" height="631" alt="Diagram depicting the DNA sequences  in the original sample as circular chromosomes, then the DNA fragmented into reads, then assembled into contigs, and then binned"/>
</a>

One way to separate contigs that belong to different species is by their taxonomic assignation. However, this can be time consuming and require a lot of computational power.
There are easier methods that perform binning to a high quality using
characteristics of the contigs, such as their GC content, their tetranucleotide frequencies (TNF), their coverage (abundance), sets of marker genes, taxonomic aligments and their preferred codons.

Most binning tools use short reads for the binning; only a few use Hi-C sequencing. Hi-C is a method of sequencing that gives spatial proximity information, as described [here](https://en.wikipedia.org/wiki/Hi-C_(genomic_analysis_technique). Different tools use different algorithms for performing the binning. A few popular tools are summarised below.  For more information see section 2.4 (Tools for metagenome binning) of [this review](https://www.sciencedirect.com/science/article/pii/S2001037021004931#s0045).

| Tool | Core algorithm | Website | Publication |
| ----------- | ----------- | -----------| -----------|
| MaxBin2      | Expectation-maximization      | http://sourceforge.net/projects/maxbin/ | [Wu et al, 2016](https://academic.oup.com/bioinformatics/article/32/4/605/1744462) |
| CONCOCT   | GaussiAN Mixture Models     | https://github.com/BinPro/CONCOCT | [Alneberg et al, 2014](https://www.nature.com/articles/nmeth.3103) |
| MetaBAT2   | Label propagation    | https://bitbucket.org/berkeleylab/metabat | [Kang et al, 2019](https://peerj.com/articles/7359/) |

There are other tools that bin MAGs using several different methods and then further refine these bins. [DAStool](https://www.nature.com/articles/s41564-018-0171-1), [MetaWRAP](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1) and Metagenome Assembled Genomes Orchestra [MAGO](https://academic.oup.com/mbe/article/37/2/593/5601623) are capable of doing this.

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) is a binning algorithm
that distinguishes between contigs that belong to different bins according to their
coverage levels and the tetranucleotide frequencies they have. We will be using this algorithm for our binning today, but first we need to prepare our assembly.

## Preparation for binning

First we index the polished reference using `bwa index`.
~~~
cd ~/cs_course/analysis/pilon/
bwa index pilon.fasta
~~~
{: .bash}

We then make a directory for the output of our binning.

~~~
cd ~/cs_course/analysis/
mkdir binning
cd binning
~~~
{: .bash}

We can then use an adapted form of the `bwa mem` command we used earlier to align our short reads to the polished assembly and determine the abundance of each contig.

~~~
( bwa mem -t 8 ../pilon/pilon.fasta ../../data/illumina_fastq/ERR2935805.fastq | samtools view - -Sb | samtools sort - -@4 -o pilon_short_read_alignment.bam ) &> binning.out &
~~~
{: .bash}
This should take around 60 minutes to complete.

Once the file is created, we can take a quick look inside.  We are already in the `binning` directory where the output file is stored so all we need to do is use the command `less` on `binning.out`.

~~~ 
less binning.out
~~~
{: .bash}

The start of the file should look like this:
~~~
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 396040 sequences (80000080 bp)...
[M::process] read 396040 sequences (80000080 bp)...
[M::mem_process_seqs] Processed 396040 reads in 97.374 CPU sec, 14.114 real sec
[M::mem_process_seqs] Processed 396040 reads in 94.057 CPU sec, 14.927 real sec
[M::process] read 396040 sequences (80000080 bp)...
[M::mem_process_seqs] Processed 396040 reads in 96.829 CPU sec, 15.646 real sec
[M::process] read 396040 sequences (80000080 bp)...
[M::mem_process_seqs] Processed 396040 reads in 92.853 CPU sec, 14.992 real sec
[M::process] read 396040 sequences (80000080 bp)...
[M::mem_process_seqs] Processed 396040 reads in 94.947 CPU sec, 15.328 real sec
[M::process] read 396040 sequences (80000080 bp)...
~~~
{: .output}

If we scroll down, the end should look like this:
~~~
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 8 ../pilon/pilon.fasta ../../../data/illumina_fastq/ERR2935805.fastq
[main] Real time: 3460.681 sec; CPU: 11813.877 sec
[bam_sort_core] merging from 12 files and 4 in-memory blocks...
~~~
{: .output}

In order to use this new BAM with MetaBAT2 we also need to index the alignment using the command `samtools index`. This only takes 2-3 minutes.

~~~
samtools index pilon_short_read_alignment.bam
~~~
{: .bash}

When we have the sorted and indexed BAM file we are then ready to use MetaBAT2.

## Binning using MetaBAT2

MetaBAT2 has been pre-installed on your instance. The [documentation](https://bitbucket.org/berkeleylab/metabat/src/master/README.md) tells us how to run the program via the command line.

The easiest way to run MetaBAT2 is using the command `runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]`. This will generate a depth file and then do the binning for us. In this example, we're also going to add the flag `-m 1500`, which sets the minimum contig length to 1500bp - any contigs shorter than this will not be binned.
~~~
runMetaBat.sh -m 1500 ../pilon/pilon.fasta pilon_short_read_alignment.bam
~~~
{: .bash}

MetaBAT2 first reads in the `.bam` file, then generates bins. This should take around 2 or 3 minutes.

While MetaBAT2 is processing the `.bam` file you will see the following output:
~~~
Executing: 'jgi_summarize_bam_contig_depths --outputDepth pilon.fasta.depth.txt --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta ../pilon/pilon.fasta pilon_short_read_alignment.bam' at Wed 28 Sep 16:45:13 BST 2022
Output depth matrix to pilon.fasta.depth.txt
Minimum percent identity for a mapped read: 0.97
minContigLength: 1000
minContigDepth: 1
Reference fasta file ../pilon/pilon.fasta
jgi_summarize_bam_contig_depths 2.15 (Bioconda) 2020-07-03T13:02:15
Output matrix to pilon.fasta.depth.txt
Reading reference fasta file: ../pilon/pilon.fasta
... 148 sequences
0: Opening bam: pilon_short_read_alignment.bam
Processing bam files
~~~
{: .output}

Once the `.bam` file has processed and binning has completed, the output will look like this:
~~~
Thread 0 finished: pilon_short_read_alignment.bam with 97425751 reads and 95337444 readsWellMapped
Creating depth matrix file: pilon.fasta.depth.txt
Closing most bam files
Closing last bam file
Finished
Finished jgi_summarize_bam_contig_depths at Wed 28 Sep 16:46:23 BST 2022
Creating depth file for metabat at Wed 28 Sep 16:46:23 BST 2022
Executing: 'metabat2  -m 1500 --inFile ../pilon/pilon.fasta --outFile pilon.fasta.metabat-bins1500-20220928_164623/bin --abdFile pilon.fasta.depth.txt' at Wed 28 Sep 16:46:23 BST 2022
MetaBAT 2 (2.15 (Bioconda)) using minContig 1500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, maxEdges 200 and minClsSize 200000. with random seed=1664379983
6 bins (14598015 bases in total) formed.
Finished metabat2 at Wed 28 Sep 16:46:25 BST 2022
~~~
{: .output}

The penultimate line tells us that MetaBAT has produced 6 bins containing 14598015 bases.

Using `ls` will show that MetaBAT2 has generated a depth file (`pilon.fasta.depth.txt`) and a directory (`pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/`). Our bins are in this directory so we should navigate into it and have a look at what files have been generated.
~~~
cd pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/
ls
~~~
{: .bash}

~~~
bin.1.fa  bin.2.fa  bin.3.fa  bin.4.fa  bin.5.fa  bin.6.fa
~~~
{: .output}

Note these output files have the file extensions of `.fa`. This is exactly the same format as a `.fasta` file but with a shortened version of the extension. See the wikipedia page on [FASTA format - file](https://en.wikipedia.org/wiki/FASTA_format#FASTA_file) for some other examples of file extensions.

Ideally we would like only one contig per bin, with a length similar the genome size of the corresponding taxa. This is challenging as this would require knowing what species are present in the mixed community, but all we have is "microbial dark matter". Instead, other statistics can demonstrate how effective our assembly and binning were.

One useful statistic is the N50 which will give an indication of the size of the contigs (fragments) each bin is made up as. In (Lesson 1)[Metagenomics 01 - 05 QC Polished Assembly](https://cloud-span.github.io/metagenomics01-qc-assembly/05-QC-polished-assembly/index.html), we got this statistic using `seqkit stats`. We can do the same again for each of the six bins.

~~~
seqkit stats -a *.fa
~~~
{: .bash}

| file     | format | type | num_seqs | sum_len | min_len | avg_len   | max_len | Q1        | Q2        | Q3        | sum_gap | N50     | Q20(%) | Q30(%) | GC(%) |
|----------|--------|------|----------|---------|---------|-----------|---------|-----------|-----------|-----------|---------|---------|--------|--------|-------|
| bin.1.fa | FASTA  | DNA  | 78       | 833410  | 3189    | 10684.7   | 28254   | 6756.0    | 8325.5    | 14036.0   | 0       | 13228   | 0.00   | 0.00   | 38.16 |
| bin.2.fa | FASTA  | DNA  | 37       | 3132462 | 4459    | 84661.1   | 334164  | 31784.0   | 59490.0   | 100708.0  | 0       | 152863  | 0.00   | 0.00   | 44.21 |
| bin.3.fa | FASTA  | DNA  | 2        | 253329  | 53561   | 126664.5  | 199768  | 53561.0   | 126664.5  | 199768.0  | 0       | 199768  | 0.00   | 0.00   | 40.99 |
| bin.4.fa | FASTA  | DNA  | 5        | 574132  | 70348   | 114826.4  | 176715  | 78563.0   | 98950.0   | 149556.0  | 0       | 149556  | 0.00   | 0.00   | 43.57 |
| bin.5.fa | FASTA  | DNA  | 2        | 6812625 | 738110  | 3406312.5 | 6074515 | 738110.0  | 3406312.5 | 6074515.0 | 0       | 6074515 | 0.00   | 0.00   | 66.18 |
| bin.6.fa | FASTA  | DNA  | 1        | 2992057 | 2992057 | 2992057.0 | 2992057 | 1496028.5 | 2992057.0 | 1496028.5 | 0       | 2992057 | 0.00   | 0.00   | 37.95 |



> ## Optional Exercise:
>
> Using the GC content and total size from the seqkit stats output, which bin might correspond to which organism in our dataset? Give your reasoning for each bin/species. (You might find it helpful to go to the section where we introduce the dataset ([Metagenomics 01 - 00 Introduction](https://cloud-span.github.io/metagenomics01-qc-assembly/00-introduction-meta/index.html)) and also do some additional reading about each species.
>
> In a later lesson we will be using a program for taxonomic assignment so you can see then how your answers compare to the results.  
>
>> ## Solution
>> - Bin 1 = maybe Enterococcus faecalis but again this one is open for interpretation as the bin is small.
>> - Bin 2 = Bacillus subtilis
>> - Bin 3 = probably Escherichia coli - but GC content is low. This one is open for interpretation as the bin is small.  
>> - Bin 4 = Saccharomyces cerevisiae (could be B subtilis but that's more abundant so probably more complete in assembly)
>> - Bin 5 = Pseudomonas aeruginosa
>> - Bin 6 = Listeria monocytogenes
> {: .solution}
{: .challenge}

> ## Recommended reading:
> Generating metagenome bins can be challenging, especially in complex community samples or where degradation of the DNA has resulted in a very incomplete assembly and short contig lengths. This workflow for binning might not work for you, and you might find that a different binning method might result in better refined MAGs for your dataset.  There are lots of other binning software methods inlcuding:
> * [CONCOCT](https://www.nature.com/articles/nmeth.3103) Clustering cOntigs with COverage and ComposiTion, the manual for running this is [here](https://github.com/BinPro/CONCOCT)
> * [Maxbin2](https://academic.oup.com/bioinformatics/article/32/4/605/1744462) uses an expectation-maximization algorithm to form bins. The link to installing maxbin2 as a conda package is [here](https://anaconda.org/bioconda/maxbin2)
> * There are also tools that can combine different binning methods and use them to refine, [DAS tool](https://www.nature.com/articles/s41564-018-0171-1) being one of them. DAS tool can also give completeness information, similarly to checkM.
> * This [review](https://www.sciencedirect.com/science/article/pii/S2001037021004931#s0055) gives a thorough look at the pros and cons of different tools used for generating MAGs and including binning .
{: .callout}

{% include links.md %}
