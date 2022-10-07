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
- "Use Metabat2 to assign the contigs to bins of different taxa."
- "Other programmes are available that are generating other bins, and these can be rationalised using tools such as DAStools"
---

## Metagenomic binning
To analyze each of the species inside our sample individually, the original genomes in the sample can be separated with a process called binning.
We call these genomes reconstructed from metagenomic assembly MAGs (Metagenome-Assembled Genomes).  
In this process, the assembled contigs from the metagenome will be assigned to different bins (FASTA files that contain certain contigs). Ideally, each bin corresponds to only one original genome (a MAG).

<a href="{{ page.root }}/fig/03-05-01.png">
  <img src="{{ page.root }}/fig/03-05-01.png" width="435" height="631" alt="Diagram depicting the DNA sequences  in the original sample as circular chromosomes, then the DNA fragmented into reads, then assembled into contigs, and then binned"/>
</a>

Although an obvious way to separate contigs that correspond to a different species is by their taxonomic assignation, this can be time consuming and require a lot of computational power.
There are easier methods that do the binning to a high quality using
characteristics of the contigs, such as their GC content, the use of tetranucleotides (composition), or their coverage (abundance).

[Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) is a binning algorithm
that distinguishes between contigs that belong to different bins according to their
coverage levels and the tetranucleotide frequencies they have.

We are going to index the polished reference first with the following command, and then use `bwa mem` command again to align our short reads to the polished assembly in order to determine the abundance of each contig.

~~~
cd analysis/pilon/
bwa index pilon.fasta
~~~
{: .bash}

We will then make a directory for the output of our binning.
~~~
cd ~/analysis/
mkdir binning
cd binning
~~~
{: .bash}

We can then run the following command (we are adapting the `bwa mem` command we've used previously).

~~~
bwa mem -t 4 ../pilon/pilon.fasta ../../data/illumina_fastq/ERR2935805.fastq | samtools view - -Sb | samtools sort - -@4 -o pilon_short_read_alignment.bam
~~~
{: .bash}

In order to use this new BAM with metabat2 we also need to index the alignment using the command `samtools index`

~~~
samtools index pilon_short_read_alignment_sort.bam
~~~
{: .bash}

When we have the sorted and indexed BAM file we are then ready to use Metabat2

Metabat2 has been pre-installed on your instance. From the documentation [README.md file](https://bitbucket.org/berkeleylab/metabat/src/master/README.md) we can see how to run Metabat2 on the command line in the section "MetaBAT 2 USAGE: running on command line".
This tells us the "easy" way to run metabat2 is using `runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]`, this will generate a depth file and then do the binning for us.
We're also going to add the flag `-m 1500`, which sets the minimum contig length it will try to bin to 1500bp
~~~
runMetaBat.sh -m 1500 ../pilon/pilon.fasta pilon_short_read_alignment.bam
~~~
{: .bash}

When you run this, `metabat2` will first read in the bam file then generate bins. This should take around 5 minutes.
You will first see the following when metabat2 is processing the bam file. It will likely stay at this point for a few minutes while the bam file is processed.
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

Once the bam file has processed and binning has completed, you should see the following output.
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

From the final line we can see that MetaBAT has produced 6 bins `6 bins (14598015 bases in total) formed.`

If you `ls` you can see that metabat2 has generated a depth file `pilon.fasta.depth.txt` and a directory `pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/`. Our bins are in this directory so we should navigate into it and have a look at what files have been generated.
~~~
cd pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/
ls
~~~
{: .bash}

~~~
bin.1.fa  bin.2.fa  bin.3.fa  bin.4.fa  bin.5.fa  bin.6.fa
~~~
{: .output}

Note these output files have the file extensions of `.fa` - this is exactly the same format as a `.fasta` file it is just a shortened version of the extension. See the wikipedia page on [FASTA format - file](https://en.wikipedia.org/wiki/FASTA_format#FASTA_file) for some other examples of file extensions.

Ideally, we would like to get only one contig per bin, with a length similar the genome size of the corresponding taxa. Since this would require knowing what species are in the mixed community along with knowing their genome size - which for many species in metagenomic samples is challenging as they make up the "microbial dark matter". As knowing this information for a lot of metagenomic samples is rare we can use other statistics that show us how good our assembly and binning have been. One useful statistic is the N50 which will give an indication of the size of the contigs (fragments) each bin is made up as. As we covered in [Metagenomics 01 - 05 QC Polished Assembly](https://cloud-span.github.io/metagenomics01-qc-assembly/05-QC-polished-assembly/index.html), we can use `seqkit stats` to get these statistics on all six of the bins.

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
> Using the GC content and total size from the seqkit stats output, which bin might correspond to which organism in our dataset? Give your reasoning for each bin / species.  
> You might find it helpful to go to the section where we introduce the dataset ([Metagenomics 01 - 00 Introduction](https://cloud-span.github.io/metagenomics01-qc-assembly/00-introduction-meta/index.html)) and also do some additional reading about each species.  
> In a later lesson we will be using a program for taxonomic assignment so you can see then how your answers compare to the results.  
>
>> ## Solution
>>   
> {: .solution}
{: .challenge}

> ## Other binning methods
> Generating metagenome bins can be challenging, especially in complex community samples or where degradation of the DNA has resulted in a very incomplete assembly and short contig lengths.  
> You might find that a different binning method might result in better refined MAGs for your dataset.  
> There are lots of other binning software methods, [CONCOCT](https://www.nature.com/articles/nmeth.3103) and [Maxbin2](https://academic.oup.com/bioinformatics/article/32/4/605/1744462) being two popular alternatives.  
> There are also tools that can combine different binning methods and use them to refine, [DAS tool](https://www.nature.com/articles/s41564-018-0171-1) being one of them. DAS tool can also give completeness information, similarly to checkM.
{: .callout}


{% include links.md %}
