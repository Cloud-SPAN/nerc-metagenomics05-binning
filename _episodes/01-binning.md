---
title: "Metagenome Binning"
teaching: 50
exercises: 10
questions:
- "How can we obtain the original genomes from a metagenome?"
objectives:
- "Obtain Metagenome-Assembled Genomes from the metagenomic assembly."
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

Although an obvious way to separate contigs that correspond to a different species is by their taxonomic assignation,
there are more reliable methods that do the binning using
characteristics of the contigs, such as their GC content, the use of tetranucleotides (composition), or their coverage (abundance).

[Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) is a binning algorithm
that distinguishes between contigs that belong to different bins according to their
coverage levels and the tetranucleotide frequencies they have.

Let's bin the sample we just assembled. The command for running Metabat2 is `runMetaBat.sh`, and the arguments it needs are the FASTA file of the assembly, the FASTQ with the forward and reverse reads, the output directory, and name. Before we run Metabat2 we need to generate a BAM file using the program BWA, and then sort these BAM files.

[BWA](http://bio-bwa.sourceforge.net/bwa.shtml) is a alignment tool, which maps reads to a reference. Due to the assembly being a *de novo* genome that we don't know what it should look like. We can instead use the assembly that we have polished, and map the reads used to generate the assembly and map them to it.

We already created a BAM file and the index in the [polishing an assembly section](https://cloud-span.github.io/metagenomics01-qc-assembly/04-polishing-assembly/index.html). However becaues there are small changes to the length of sequence of the reference after using pilon, we will need to generate a new BAM fileand index.

We are going to index the polished reference first with the following command, and then use bwa mem command again.
~~~
bwa index pilon.fasta

bwa mem -t 4  pilon.fasta ../../ERR3152367_sub5_filtered.fastq | samtools view - -Sb | samtools sort - -@4 -o pilon_short_read_alignment.bam >> pilon_alignment.out 2>&1  &
~~~
{: .bash}

In order to use this new BAM with metabat2 we also need to sort the order of the alignments using the command `samtools sort`.


We sort the BAM file we generated in the last lesson
~~~
samtools sort -o pilon_short_read_alignment_sort.bam pilon_short_read_alignment.bam

We then index this BAM file
samtools index pilon_short_read_alignment_sort.bam
~~~
{: .bash}

When we have the sorted BAM file we are then ready to use Metabat2

~~~
runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
runMetaBat.sh assembly.fasta sample1.bam &
~~~
{: .bash}


~~~
 cd ~/analysis/
 mkdir Metabat2
 cd Metabat2
 runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
 runMetaBat.sh pilon.fasta pilon_short_read_alignment_sort.bam &
~~~
{: .bash}

Once you have ran metabat2 you should have this output to the screen

~~~
Executing: 'jgi_summarize_bam_contig_depths --outputDepth pilon.fasta.depth.txt --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta pilon.fasta pilon_short_read_alignment_sort.bam' at Wed  7 Sep 17:54:51 BST 2022
Output depth matrix to pilon.fasta.depth.txt
Minimum percent identity for a mapped read: 0.97
minContigLength: 1000
minContigDepth: 1
Reference fasta file pilon.fasta
jgi_summarize_bam_contig_depths 2.15 (Bioconda) 2020-01-04T21:10:40
Output matrix to pilon.fasta.depth.txt
Reading reference fasta file: pilon.fasta
... 146 sequences
0: Opening bam: pilon_short_read_alignment_sort.bam
Processing bam files
Thread 0 finished: pilon_short_read_alignment_sort.bam with 928919 reads and 53364 readsWellMapped
Creating depth matrix file: pilon.fasta.depth.txt
Closing most bam files
Closing last bam file
Finished
Finished jgi_summarize_bam_contig_depths at Wed  7 Sep 17:55:17 BST 2022
Creating depth file for metabat at Wed  7 Sep 17:55:17 BST 2022
Executing: 'metabat2  --inFile pilon.fasta --outFile pilon.fasta.metabat-bins-20220907_175517/bin --abdFile pilon.fasta.depth.txt' at Wed  7 Sep 17:55:17 BST 2022
MetaBAT 2 (2.15 (Bioconda)) using minContig 2500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, maxEdges 200 and minClsSize 200000. with random seed=1662569717
0 bins (0 bases in total) formed.
Finished metabat2 at Wed  7 Sep 17:55:20 BST 2022
~~~
{: .output}

If you would like to see other options to alter your output you can get these by accepting the metabat2 manual

~~~
metabat2 -h

Allowed options:
  -h [ --help ]                     produce help message
  -i [ --inFile ] arg               Contigs in (gzipped) fasta file format [Mandatory]
  -o [ --outFile ] arg              Base file name and path for each bin. The default output is fasta format.
                                    Use -l option to output only contig names [Mandatory].
  -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited;
                                    the first column should be contig names, and the first row will be
                                    considered as the header and be skipped) [Optional].
  -m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
  --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
                                    among contigs. The greater, the more sensitive.
  --minS arg (=60)                  Minimum score of a edge for binning (should be between 1 and 99). The
                                    greater, the more specific.
  --maxEdges arg (=200)             Maximum number of edges per node. The greater, the more sensitive.
  --pTNF arg (=0)                   TNF probability cutoff for building TNF graph. Use it to skip the
                                    preparation step. (0: auto).
  --noAdd                           Turning off additional binning for lost or small contigs.
  --cvExt                           When a coverage file without variance (from third party tools) is used
                                    instead of abdFile from jgi_summarize_bam_contig_depths.
  -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig in each library for binning.
  --minCVSum arg (=1)               Minimum total effective mean coverage of a contig (sum of depth over
                                    minCV) for binning.
  -s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
  -t [ --numThreads ] arg (=0)      Number of threads to use (0: use all cores).
  -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
  --saveCls                         Save cluster memberships as a matrix format
  --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
  --noBinOut                        No bin output. Usually combined with --saveCls to check only contig
                                    memberships
  --seed arg (=0)                   For exact reproducibility. (0: use random seed)
  -d [ --debug ]                    Debug output
  -v [ --verbose ]                  Verbose output
  ~~~
{: .output}

In MetaBAT 2, parameter optimization will be unnecessary, though we allowed a few parameters so that advanced users might play with them.

You can decrease -m (--minContig) when the qualities of both assembly and formed bins with default value are very good.

'-i' input file should be either fasta or gzipped fasta file. (since v0.32.3)

-p option is for utilizing paired info from short reads. It may improve sensitivity.
'--p1' and '--p2' should be both high to maintain great specificity. Usually p1 >= p2 performs better.
--minProb mainly controls the scope and sensitivity of binning. A smaller number improves sensitivity. It should be < p1, p2.

--minBinned mainly controls the specificity. A greater number improves specificity. Usually <= 50.
Use --verysensitive on simple community for more inclusive binning.

--minCorr would include contigs which are closely correlated in abundance but somewhat different in absolute abundance. More effective in availability of many samples.
Recruiting by correlation would be activated only if # of samples >= minSamples and be disabled (for better specificity) when minContigByCorr > minContig.

Smaller contigs (>1000) will be given a chance to be recruited to existing bins when # of samples >= minSamples by default setting.

'--saveDistance' option saves a lot of computations when multiple binning attempts are executed with different parameter settings.

'--unbinned' option generates a file for unbinned contigs.

'-B' option is for ensemble binning. Recommended to be 20 or more. Should be >= 10 for reasonable results. It tends to generate reduced number of better quality bins at the cost of some additional computation.

'--pB' option controls for sensitivity and specificity tradeoff in ensemble binning. The smaller, the sensitive. Range is between 0 to 100. The default is 50.

Produced bins would be stochastic if ensemble binning was used. --seed would minimize the stochasticity but still there would be slight difference.

Each bin will be saved as a fasta format


{% include links.md %}
