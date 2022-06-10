---
title: "Metagenome Binning"
teaching: 50
exercises: 10
questions:
- "How can we obtain the original genomes from a metagenome?"
objectives:
- "Obtain Metagenome-Assembled Genomes from the metagenomic assembly."
- "Check the quality of the Metagenome-Assembled genomes."  
keypoints:
- "Metagenome-Assembled Genomes (MAGs) sometimes are obtained from curated contigs grouped into bins."
- "Use Metabat2 to assign the contigs to bins of different taxa."
- "Use CheckM to evaluate the quality of each Metagenomics-Assembled Genome."
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

Let's bin the sample we just assembled. The command for running MaxBin is `runMetaBat.sh`, and the arguments it needs are the FASTA file of the assembly, the FASTQ with the forward and reverse reads, the output directory, and name. Before we run Metabat2 we need to generate a BAM file using the program BWA, and then sort these BAM files.

[BWA](http://bio-bwa.sourceforge.net/bwa.shtml) is a alignment tool, which maps reads to a reference. Due to the assembly being a *de novo* genome that we don't know what it should look like. We can instead use the assembly that we have polished, and map the reads used to generate the assembly and map them to it.

~~~
$ # We use bwa index to generate an inde for the fasta containing the contigs
$ bwa index assembly.fa
$ # We then align the reads against contigs and sort this BAM file
$ bwa mem  -t 16 assembly.fa shortreads.fq.gz | samtools sort -o alignment.bam
$ # We then index this BAM file
$ samtools index alignment.bam
~~~
{: .bash}

When we have the sorted BAM file we are then ready to use Metabat2

~~~
$ runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
$ runMetaBat.sh assembly.fasta sample1.bam &
~~~
{: .bash}


~~~
$ cd ~/dc_workshop/results/assembly_JC1A
$ mkdir Metabat2
$ runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
$ runMetaBat.sh assembly.fasta sample1.bam &
~~~
{: .bash}
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


## Quality check

The quality of a MAG is highly dependent on the size of the genome of the species, its abundance
in the community, and the depth at which we sequenced it.
Two important things that can be measured to know its quality are completeness (is the MAG a complete genome?)
and if it is contaminated (does the MAG contain only one genome?).

[CheckM](https://github.com/Ecogenomics/CheckM) is a good program to see the quality of our MAGs.
It gives a measure of the completeness and the contamination by counting marker genes in the MAGs.
The lineage workflow that is a part of CheckM places your bins in a reference tree to know to which lineage it corresponds to and to use the appropriate marker genes to estimate the quality parameters. Unfortunately, the lineage workflow uses a lot of memory so it can't run in our machines, but we can tell CheckM to use marker genes from Bacteria only, to spend less memory.
This is a less accurate approach but it can also be very useful if you want all of your bins analyzed with the same markers.

We will run the taxonomy workflow specifying the use of markers at the domain level, specific for the rank Bacteria,
we will specify that our bins are in FASTA format, that they are located in the `MAXBIN` directory
and that we want our output in the `CHECKM/` directory.
~~~
$ mkdir CHECKM
$ checkm taxonomy_wf domain Bacteria -x fasta MAXBIN/ CHECKM/
~~~
{: .bash}

The run will end with our results printed in the console.
~~~
--------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id     Marker lineage   # genomes   # markers   # marker sets   0    1    2    3    4   5+   Completeness   Contamination   Strain heterogeneity  
--------------------------------------------------------------------------------------------------------------------------------------------------------
  JP4D.002      Bacteria         5449        104            58        3    34   40   21   5   1       94.83           76.99              11.19          
  JP4D.004      Bacteria         5449        104            58        12   40   46   6    0   0       87.30           51.64               3.12          
  JP4D.001      Bacteria         5449        104            58        24   65   11   3    1   0       70.48           13.09               0.00          
  JP4D.003      Bacteria         5449        104            58        44   49   11   0    0   0       64.44           10.27               0.00          
--------------------------------------------------------------------------------------------------------------------------------------------------------

~~~
{: .output}

To have these values in an output that is more usable and shearable we can now run the quality step of CheckM `checkm qa`
and make it print the output in a `TSV` table, instead of the console. In this step, we can ask CheckM to give us more parameters, like contig number and length.

Ideally, we would like to get only one contig per bin, with a length similar the genome size of the corresponding taxa. Since this scenario is very difficult to obtain we can use parameters that show us how good is our assembly. Here are some of the most common metrics:
If we arrange our contigs by size, from larger to smaller, and divide the whole sequence in half, N50 is the size of the smallest contig in the half that has the larger contigs; and L50 is the number of contigs in this half of the sequence. So we want big N50 and small L50 values for our genomes. Read [What’s N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/).

To get the table with these extra parameters we need to specify the file of the markers that CheckM used in the previous step `Bacteria.ms`, the name of the output file we want `quality_JP4D.tsv`, that we want a table `--tab_table`, and the option number 2 `-o 2` is to ask for the extra parameters printed on the table.
~~~
$  checkm qa CHECKM/Bacteria.ms CHECKM/ --file CHECKM/quality_JP4D.tsv --tab_table -o 2
~~~
{: .bash}
The table we just made looks like [this](https://github.com/carpentries-incubator/metagenomics/blob/gh-pages/files/quality_JP4D.tsv).
This will be very useful when you need to document your work or communicate it.

The question of, how much contamination we can tolerate and how much completeness do we need, certainly depends on the scientific question being tackled, but in the [CheckM](https://genome.cshlp.org/content/25/7/1043) paper, there are some parameters that we can follow.


> ## Exercise 1: Discuss the quality of the obtained MAGs
>
> Fill in the blanks to complete the code you need to download the `quality_JP4D.tsv` to your local computer:
> ~~~
> ____ dcuser____ec2-18-207-132-236.compute-1.amazonaws.com____/home/dcuser/dc_workshop/mags/CHECKM/quality_JP4D.tsv ____
> ~~~
> {: .bash}
>
>> ## Solution
>>In a terminal that is standing on your local computer do:
>> ```
>>$ scp dcuser@ec2-18-207-132-236.compute-1.amazonaws.com:/home/dcuser/dc_workshop/mags/CHECKM/quality_JP4D.tsv <the destination directory of your choice>
>> ```
>>{: .bash}
>>
> {: .solution}
> Then open the table in a spreadsheet and discuss with your team which of the parameters in the table do you find useful.
{: .challenge}

{: .bash}
{% include links.md %}
