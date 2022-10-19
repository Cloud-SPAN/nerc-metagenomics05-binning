---
title: "QC of metagenome bins"
teaching: 50
exercises: 10
questions:
- "How can we assess the quality of the metagenome bins?"
objectives:
- "Check the quality of the Metagenome-Assembled Genomes (MAGs)."
- "Understanding MIMAG quality standards."  
keypoints:
- "CheckM can be used to evaluate the quality of each Metagenomics-Assembled Genome."
- "We can use the percentage contamination and completion to identify the quality of these bins."
- "There are MIMAG standards which can be used to categorise the quality of a MAG."
- "Many MAGs will be incomplete, but that does not mean that this data is not still useful for downstream analysis."
---

## Quality check

The quality of a metagenome-assembled genome (MAG) or bin is highly dependent on several things:
- the depth of sequencing
- the abundance of the organism in the community
- how successful the assembly was
- how successful the polishing (if used) was

In order to determine the quality of a MAG we can look at two different metrics. These are:
1. completeness (i.e. how much of the genome is captured in the MAG?) and 
2. contamination (i.e. do all the sequences in the MAG belong to the same organism?).

We can use the program [CheckM](https://github.com/Ecogenomics/CheckM) to determine the quality of MAGs. CheckM uses a collection of domain and lineage-specific markers to estimate completeness and contamination of a MAG. This [short YouTube video](https://youtu.be/sLtSDs3sh6k) by Dr Robert Edwards explains how CheckM uses a hidden Markov model to calculate the level of contamination and completeness of bins, based on marker gene sets.

CheckM has multiple different workflows available which are appropriate for different datasets. See [CheckM documentation on Workflows](https://github.com/Ecogenomics/CheckM/wiki/Workflows) for more information.

We will be using the lineage-specific workflow here. `lineage_wf` places your bins in a reference tree to determine which lineage it corresponds to. This allows it to use the appropriate marker genes to estimate quality parameters.

First let's move into our analysis folder.
~~~
cd ~/cs_course/analysis/
~~~
{: .bash}

CheckM has been pre-installed on the instance so we can check the help documentation for the lineage-specific workflow using the `'h` tag..  

~~~
checkm lineage_wf -h
~~~
{: .bash}

> ## CheckM help documentation
> ~~~
> usage: checkm lineage_wf [-h] [-r] [--ali] [--nt] [-g] [-u UNIQUE] [-m MULTI]
>                          [--force_domain] [--no_refinement]
>                          [--individual_markers] [--skip_adj_correction]
>                          [--skip_pseudogene_correction]
>                          [--aai_strain AAI_STRAIN] [-a ALIGNMENT_FILE]
>                          [--ignore_thresholds] [-e E_VALUE] [-l LENGTH]
>                          [-f FILE] [--tab_table] [-x EXTENSION] [-t THREADS]
>                          [--pplacer_threads PPLACER_THREADS] [-q]
>                          [--tmpdir TMPDIR]
>                          bin_input output_dir
>
> Runs tree, lineage_set, analyze, qa
>
> positional arguments:
>   bin_input             directory containing bins (fasta format) or path to file describing genomes/genes - tab separated in 2 or 3 columns [genome ID, genome fna, genome translation file (pep)]
>   output_dir            directory to write output files
>
> optional arguments:
>   -h, --help            show this help message and exit
>   -r, --reduced_tree    use reduced tree (requires <16GB of memory) for determining lineage of each bin
>   --ali                 generate HMMER alignment file for each bin
>   --nt                  generate nucleotide gene sequences for each bin
>   -g, --genes           bins contain genes as amino acids instead of nucleotide contigs
>   -u, --unique UNIQUE   minimum number of unique phylogenetic markers required to use lineage-specific marker set (default: 10)
>   -m, --multi MULTI     maximum number of multi-copy phylogenetic markers before defaulting to domain-level marker set (default: 10)
>   --force_domain        use domain-level sets for all bins
>   --no_refinement       do not perform lineage-specific marker set refinement
>   --individual_markers  treat marker as independent (i.e., ignore co-located set structure)
>   --skip_adj_correction
>                         do not exclude adjacent marker genes when estimating contamination
>   --skip_pseudogene_correction
>                         skip identification and filtering of pseudogenes
>   --aai_strain AAI_STRAIN
>                         AAI threshold used to identify strain heterogeneity (default: 0.9)
>   -a, --alignment_file ALIGNMENT_FILE
>                         produce file showing alignment of multi-copy genes and their AAI identity
>   --ignore_thresholds   ignore model-specific score thresholds
>   -e, --e_value E_VALUE
>                         e-value cut off (default: 1e-10)
>   -l, --length LENGTH   percent overlap between target and query (default: 0.7)
>   -f, --file FILE       print results to file (default: stdout)
>   --tab_table           print tab-separated values table
>   -x, --extension EXTENSION
>                         extension of bins (other files in directory are ignored) (default: fna)
>   -t, --threads THREADS
>                         number of threads (default: 1)
>   --pplacer_threads PPLACER_THREADS
>                         number of threads used by pplacer (memory usage increases linearly with additional threads) (default: 1)
>   -q, --quiet           suppress console output
>   --tmpdir TMPDIR       specify an alternative directory for temporary files
>
> Example: checkm lineage_wf ./bins ./output
> ~~~
> {: .output}
{: .solution}

This readout tells us what we need to include in the command:
- the directory that contains the bins (`pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/)` 
- the `x` flag telling CheckM the format of our bins
- the `--reduced_tree` flag to limit the memory requirements 
- the `-f` flag to specify an output file name/format
- the `--tab_table` flag  so the output is in a tab-separated format
- the  `-t` flag to set the number of threads used to four, which is the number we have on our instance

As a result our command looks like this:
~~~
checkm lineage_wf  -x fa mbinning/pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/ checkm/ --reduced_tree -t 8 --tab_table -f MAGs_checkm.tsv
~~~
{: .bash}

When the run ends our results we can open our results file.
~~~
less MAGs_checkm.tsv
~~~
{: .bash}

~~~
Bin Id  Marker lineage  # genomes       # markers       # marker sets   0       1       2       3       4       5+      Completeness    Contamination   Strain heterogeneity
bin.1   k__Bacteria (UID203)    5449    104     58      100     4       0       0       0       0       2.19    0.00    0.00
bin.2   k__Bacteria (UID203)    5449    104     58      57      47      0       0       0       0       70.69   0.00    0.00
bin.3   root (UID1)     5656    56      24      56      0       0       0       0       0       0.00    0.00    0.00
bin.4   g__Bacillus (UID864)    93      711     241     595     116     0       0       0       0       6.42    0.00    0.00
bin.5   o__Pseudomonadales (UID4488)    185     813     308     1       807     5       0       0       0       99.68   0.61    0.00
bin.6   c__Bacilli (UID285)     586     325     181     1       324     0       0       0       0       99.45   0.00    0.00
~~~
{: .output}

Running this workflow is equivalent to running six separate CheckM commands. The [CheckM documentation](https://github.com/Ecogenomics/CheckM/wiki/Workflows) explains this is more detail.

> ## Exercise 1: Downloading the tsv file.
>
> Fill in the blanks to complete the code you need to download the `MAGs_checkm.tsv` to your local computer using SCP:
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
{: .challenge}

How much contamination we can tolerate and how much completeness we need depends on the scientific question being tackled.

To help us, we can use a standard called Minimum Information about a Metagenome-Assembled Genome (MIMAG), developed by the Genomics Standard Consortium. You can read more about MIMAG in [this 2017 paper](https://www.nature.com/articles/nbt.3893).

As part of the standard,a framework to determine MAG quality from statistics is outlined. A MAG can be assigned one of three different metrics: High, Medium or Low quality draft metagenome assembled genomes.

See the table below for an overview of each category.

| Quality Category | Completeness | Contamination | rRNA/tRNA encoded|
| ----------- | ----------- | -----------| -----------|
| High      | > 90%       | ≤ 5% | Yes (≥ 18 tRNA and all rRNA)|
| Medium   | ≥ 50%        | ≤ 10% | No |
| Low   | < 50%      | ≤ 10%| No |

We have already determined the **completeness** and **contamination** of each of our MAGs using CheckM. Next we will use a program to determine which rRNA and tRNAs are present in each MAG. 

Note that due to the difficulty in assembly of short-read metagenomes, often just a completeness of >90% and a contamination of ≤ 5% is treated as a good quality MAG.

> ## Exercise 2: Explore the quality of the obtained MAGs
>
> Once you have downloaded the `MAGs_checkm.tsv` file, you can open it in Excel or another spreadsheet program. If you didn't manage to download the file, or do not have an appropriate program to view it in you can see or download our example file [here](https://drive.google.com/file/d/11vJr5uHBx6J57sazLoU7oOl1peB9OxxJ/view?usp=sharing).
>
> Looking at the results of our quality checks, what category would each of our MAGs fall into (ignore the tRNA and rRNA requirement for now)?
> {: .bash}
>
>> ## Solution
>> There are two potential high quality draft metagenome assembled genomes in bin.5 and bin.6. We also have one medium quality draft MAG in bin.2. Finally, there are three low quality MAGs in bin.1, bin.3 and bin.4.
>>
>> Your bins may have different names/numbers to these but you should still see similar results.
>>
> {: .solution}
{: .challenge}

{% include links.md %}
