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

~~~
cd ~/analysis/
~~~
{: .bash}

CheckM has been pre-installed on the instance so we can check the help documentation for the lineage_wf.  

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

From this we can see that along with the directory that contains the bins and an output directory we also need to add a flag to tell CheckM the format of our bins, which in our case is `.fa`. We also need to pass the `--reduced_tree` flag which limits the memory requirements so we can successfully run it on the amount of compute we have. The flag `-f` will specify and output file and we also want to use the flag `--tab_table` so the output is in a tab-separated format. Finally, we want to use `-t 4` to set the number of threads used to four, which is the number we have on our instance.

~~~
checkm lineage_wf  -x fasta metabat2/pilon.fasta.metabat-bins1500-YYYMMDD_HHMMSS/ checkm/ --reduced_tree -t 4 --tab_table -f MAGs_checkm.tsv
~~~
{: .bash}

The run will end with our results in a file called which we can open with `less MAGs_checkm.tsv`
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

Running this workflow means that we run four checkM commands in one go rather than individually. This is outlined [here](https://github.com/Ecogenomics/CheckM/wiki/Workflows), which shows lineage_wf is equivalent to running the following commands:

~~~
  checkm tree <bin folder> <output folder>
  checkm tree_qa <output folder>
  checkm lineage_set <output folder> <marker file>
  checkm analyze <marker file> <bin folder> <output folder>
  checkm qa <marker file> <output folder>
  checkm qa checkM/checkM.ms checkM/ --file checkM/MAGs_checkm.tsv --tab_table -o 2
~~~
{: .output}

> ## Exercise 1: Downloading the tsv file.
>
> Fill in the blanks to complete the code you need to download the `MAGs_checkm.tsv` to your local computer:
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

The question of how much contamination we can tolerate and how much completeness do we need depends a lot on the scientific question being tackled.

In the Minimum information about a metagenome-assembled genome (MIMAG) microbial standards [Bowers, R., Kyrpides, N., Stepanauskas, R. et al., 2017](https://www.nature.com/articles/nbt.3893). A framework to determine MAG quality from statistics is outlined. Three different metrics to be assigned as either; High, Medium or Low quality draft metagenome assembled genomes.
See the table below for an overview of each category.

| Quality Category | Completeness | Contamination | rRNA/tRNA encoded|
| ----------- | ----------- | -----------| -----------|
| High      | > 90%       | ≤ 5% | Yes (≥ 18 tRNA and all rRNA)|
| Medium   | ≥ 50%        | ≤ 10% | No |
| Low   | < 50%      | ≤ 10%| No |

Using CheckM we have determined the Completeness and Contamination of each of our MAGs. We will be using a program in the next episode to determine which rRNA and tRNAs are present in each MAG. However, due to the difficulty in assembly of short-read metagenomes often just a completeness >90% and a contamination ≤ 5% is treated as a good quality MAG.

> ## Exercise 2: Explore the quality of the obtained MAGs
>
> Once you have downloaded the `MAGs_checkm.tsv` file, you can open it in Excel or another spreadsheet program. If you didn't manage to download the file, or do not have an appropriate program to view it in you can see or download our example file [here](https://drive.google.com/file/d/11vJr5uHBx6J57sazLoU7oOl1peB9OxxJ/view?usp=sharing).
> Looking at this file what category would each of our MAGs fall into (we will ignore the tRNA and rRNA requirement for now)
> {: .bash}
>
>> ## Solution
>> From the file above, we can see there are two potential high quality draft metagenome assembled genomes which are bin.5 and bin.6. We also have one medium quality draft MAG in bin.2. Finally there are three low quality MAGs in bin.1, bin.3 and bin.4.
Your bins may have a different name/number to these but you should still have seen similar results. 
>>
> {: .solution}
{: .challenge}

{% include links.md %}
