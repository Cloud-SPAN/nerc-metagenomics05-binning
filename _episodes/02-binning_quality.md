---
title: "QC of metagenome bins"
teaching: 50
exercises: 10
questions:
- "How can be assess the quality of the metagenome bins?"
objectives:
- "Check the quality of the Metagenome-Assembled genomes."
- "Understanding MIMAG quality standards."  
keypoints:
- "CheckM can be used to evaluate the quality of each Metagenomics-Assembled Genome."
- "We can use the percentage contamination and completion to identify the quality of these bins."
- "There are MIMAG standards which can be used to categorise the quality of a MAG."
- "Many MAGs will be incomplete, but that does not mean that this data is not still useful for downstream analysis."
---

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
we will specify that our bins are in FASTA format, that they are located in the `Metabat2` directory
and that we want our output in the `CHECKM/` directory.
~~~
$ mkdir CHECKM
$ checkm taxonomy_wf domain Bacteria -x fasta Metabat2/ CHECKM/
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
