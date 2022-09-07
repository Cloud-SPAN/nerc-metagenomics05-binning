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
in the community, and the depth at which we sequenced it.Two important things that can be measured to know its quality are completeness (is the MAG a complete genome?) and if it is contaminated (does the MAG contain only one genome?).

[CheckM](https://github.com/Ecogenomics/CheckM) is a good program to see the quality of our MAGs.
It gives a measure of the completeness and the contamination by counting marker genes in the MAGs. Here is a [short youtube video](https://youtu.be/sLtSDs3sh6k) by Dr Robert Edwards that explains how CheckM uses a hidden Markov model to calculate the level of contamination and completeness of bins, based on marker gene sets.

Two of the main methods are the lineage workflow `lineage_wf` and the taxnonmic workflow `taxonomy_wf`. The lineage workflow places your bins in a reference tree to know to which lineage it corresponds to and to use the appropriate marker genes to estimate the quality parameters. The taxonomy workflow works similarly but can be used for generating more specific markers for a specific taxonomic group such as a specific phylum. If you are concerned which is more appropriate for you the [checkM](https://github.com/Ecogenomics/CheckM) manual compared results across 1000 randomly selected genomes using [GTDB representative genomes](https://gtdb.ecogenomic.org/), and identical results were obtained using the lineage_wf, taxonomy_wf and ssu_finder methods. Details on how to run both workflows are given [here](https://github.com/Ecogenomics/CheckM/wiki/Workflows).



We will run the lineage workflow and will specify that our bins are in FASTA format, that they are located in the `Metabat2` directory and that we want our output in the `checkM/` directory. We will be using the `reduced_tree` option to keep our RAM consumption to 16Gb, and `-t 4` to set the number of thread to 4 because these are available on the instance and will speed up the process. If you are on a High performance computing cluster (HPC), you can remove the reduced_tree option.
~~~
$ mkdir checkM
$ checkm lineage_wf  -x fasta Metabat2/ checkM/ --reduced_tree -t 4 -o 2 --tab_table -f MAGs_checkm.tsv
~~~
{: .bash}

The run will end with our results in a file called which we can open with `less MAGs_checkm.tsv`
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

Running this workflow means that we run four checkM commands in one go rather than individually. This is outlined [here](https://github.com/Ecogenomics/CheckM/wiki/Workflows), which shows lineage_wf is equivalent to running the following commands:

~~~
  checkm tree <bin folder> <output folder>
  checkm tree_qa <output folder>
  checkm lineage_set <output folder> <marker file>
  checkm analyze <marker file> <bin folder> <output folder>
  checkm qa <marker file> <output folder>
  checkm qa checkM/checkM.ms checkM/ --file checkM/MAGs_checkm.tsv --tab_table -o 2
~~~
{: .bash}


Ideally, we would like to get only one contig per bin, with a length similar the genome size of the corresponding taxa. Since this scenario is very difficult to obtain we can use parameters that show us how good is our assembly. Here are some of the most common metrics:
If we arrange our contigs by size, from larger to smaller, and divide the whole sequence in half, N50 is the size of the smallest contig in the half that has the larger contigs; and L50 is the number of contigs in this half of the sequence. So we want big N50 and small L50 values for our genomes. Read [What’s N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/).

To get the table with these extra parameters we need to specify the file of the markers that CheckM used in the previous step `checkM.ms`, the name of the output file we want `MAGs_checkm.tsv`, that we want a table `--tab_table`, and the option number 2 `-o 2` is to ask for the extra parameters printed on the table. This means we will run the checkm qa part of the workflow again seperately with these addition commands:

~~~
  checkm qa checkM/checkM.ms checkM/ --file checkM/MAGs_checkm.tsv --tab_table -o 2
~~~
{: .bash}

The question of, how much contamination we can tolerate and how much completeness do we need, certainly depends on the scientific question being tackled, but in the [CheckM](https://genome.cshlp.org/content/25/7/1043) paper, there are some parameters that we can follow.

## Other binning methods

Generating metagenome bins can be challenging, especially in complex community samples or where degradation of the DNA has resulted in a very incomplete assembly and short contig lengths. You might find that a different binning method might result in better refined MAGs for your dataset. There are lots of other binning software methods, [CONCOCT](https://www.nature.com/articles/nmeth.3103) and [Maxbin2](https://academic.oup.com/bioinformatics/article/32/4/605/1744462) being two popular alternatives. There are also tools that can combine different binning methods and use them to refine, [DAS tool](https://www.nature.com/articles/s41564-018-0171-1) being one of them. DAS tool can also give completeness information, similarly to checkM.

> ## Exercise 1: Explore the quality of the obtained MAGs
>
> Fill in the blanks to complete the code you need to download the `quality_JP4D.tsv` to your local computer:
> ~~~
> ____ csuser____ec2-18-207-132-236.compute-1.amazonaws.com____/home/csuser/cs_workshop/mags/checkM/MAGs_checkm.tsv ____
> ~~~
> {: .bash}
>
>> ## Solution
>>In a terminal logged into your local machine type:
>> ```
>>$ scp csuser@ec2-18-207-132-236.compute-1.amazonaws.com:/home/csuser/cs_workshop/mags/checkM/MAGs_checkm.tsv <the destination directory of your choice>
>> ```
>>{: .bash}
>>
> {: .solution}
> Open the spreadsheet in excel and think about  which of the parameters in the table you find useful.
{: .challenge}

{: .bash}
{% include links.md %}
