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

The quality of a metagenome-assembled genome (MAG) or bin is highly dependent on the depth of sequencing, the abundance of the organism in the community and and how successful the assembly and polishing (if any) was.

In order to determine the quality of a MAG we can look at two different metrics. These are the completeness (i.e. how much of the genome is captured in the MAG) and contamination (i.e. do all the sequences in the MAG belong to the same organism).

We can use the program [CheckM](https://github.com/Ecogenomics/CheckM) to determine the quality of MAGs. CheckM uses a collection of domain and lineage-specific markers to estimate completeness and contamination of a MAG. Here is a [short youtube video](https://youtu.be/sLtSDs3sh6k) by Dr Robert Edwards that explains how CheckM uses a hidden Markov model to calculate the level of contamination and completeness of bins, based on marker gene sets.

CheckM has multiple different workflows available which are appropriate for different datasets, see [CheckM documentation on Workflows](https://github.com/Ecogenomics/CheckM/wiki/Workflows) for more information.

We will be using the lineage workflow here. `lineage_wf` places your bins in a reference tree to determine which lineage it corresponds to in order to use the appropriate marker genes to estimate the quality parameters.

CheckM has been pre-installed on the instance so we can first check the help documentation.  

> ## CheckM help documentation
> ~~~
>
>                 ...::: CheckM v1.2.1 :::...
>
>   Lineage-specific marker set:
>     tree         -> Place bins in the reference genome tree
>     tree_qa      -> Assess phylogenetic markers found in each bin
>     lineage_set  -> Infer lineage-specific marker sets for each bin
>
>   Taxonomic-specific marker set:
>     taxon_list   -> List available taxonomic-specific marker sets
>     taxon_set    -> Generate taxonomic-specific marker set
>
>   Apply marker set to genome bins:
>     analyze      -> Identify marker genes in bins
>     qa           -> Assess bins for contamination and completeness
>
>   Common workflows (combines above commands):
>     lineage_wf   -> Runs tree, lineage_set, analyze, qa
>     taxonomy_wf  -> Runs taxon_set, analyze, qa
>
>   Reference distribution plots:
>     gc_plot      -> Create GC histogram and delta-GC plot
>     coding_plot  -> Create coding density (CD) histogram and delta-CD plot
>     tetra_plot   -> Create tetranucleotide distance (TD) histogram and delta-TD plot
>     dist_plot    -> Create image with GC, CD, and TD distribution plots together
>
>   General plots:
>     nx_plot      -> Create Nx-plots
>     len_hist     -> Sequence length histogram
>     marker_plot  -> Plot position of marker genes on sequences
>     gc_bias_plot -> Plot bin coverage as a function of GC
>
>   Bin exploration and modification:
>     unique       -> Ensure no sequences are assigned to multiple bins
>     merge        -> Identify bins with complementary sets of marker genes
>     outliers     -> [Experimental] Identify outlier in bins relative to reference distributions
>     modify       -> [Experimental] Modify sequences in a bin
>
>   Utility functions:
>     unbinned     -> Identify unbinned sequences
>     coverage     -> Calculate coverage of sequences
>     tetra        -> Calculate tetranucleotide signature of sequences
>     profile      -> Calculate percentage of reads mapped to each bin
>     ssu_finder   -> Identify SSU (16S/18S) rRNAs in sequences
>
>   Use 'checkm data setRoot <checkm_data_dir>' to specify the location of CheckM database files.
>
>   Usage: checkm <command> -h for command specific help
>
> ~~~
> {: .output}
{: .solution}

We will run the lineage workflow and will specify that our bins are in FASTA format, that they are located in the `Metabat2` directory and that we want our output in the `checkM/` directory. We will be using the `reduced_tree` option to keep our RAM consumption to 16Gb, and `-t 4` to set the number of thread to 4 because these are available on the instance and will speed up the process. If you are on a High performance computing cluster (HPC), you can remove the reduced_tree option.
~~~
cd ~/analysis/
mkdir checkM
checkm lineage_wf  -x fasta Metabat2/ checkM/ --reduced_tree -t 4 -o 2 --tab_table -f MAGs_checkm.tsv
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
