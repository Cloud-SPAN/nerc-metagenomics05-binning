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

Prokka produces multiple different file types which you can see in the table below. We are mainly interested in `.faa` and `.tsv` but many of the other files are useful for submission to different databases.

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
| .tsv	 | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product |

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

If we navigate into the output file created by Prokka we can then list and see that Prokka has generated many files.

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

As we said earlier the two files we are mainly interested in are those with the extension `.tsv` and `.faa`. The tsv file will contain information about every gene identified by Prokka, including its length and name. The faa file is a FASTA file containing the amino acid sequence of every gene that has been identified.

We can take a look at the output of both using `head`
~~~
head bin.6.tsv
~~~
{: .bash}

~~~
locus_tag	ftype	length_bp	gene	EC_number	COG	product
JODGFBLK_00001	CDS	768	lysS	6.1.1.6		Lysine--tRNA ligase
JODGFBLK_00002	CDS	1002	dusB	1.3.1.-	COG0042	tRNA-dihydrouridine synthase B
JODGFBLK_00003	CDS	219				hypothetical protein
JODGFBLK_00004	CDS	504	folK	2.7.6.3	COG0801	2-amino-4-hydroxy-6-hydroxymethyldihydropteridine pyrophosphokinase
JODGFBLK_00005	CDS	363	folB	4.1.2.25	COG1539	Dihydroneopterin aldolase
JODGFBLK_00006	CDS	858	folP	2.5.1.15		Dihydropteroate synthase
JODGFBLK_00007	CDS	882	pabC	4.1.3.38	COG0115	Aminodeoxychorismate lyase
JODGFBLK_00008	CDS	585	pabA	2.6.1.85	COG0512	Aminodeoxychorismate/anthranilate synthase component 2
JODGFBLK_00009	CDS	1410	pabB	2.6.1.85	COG0147	Aminodeoxychorismate synthase component 1
~~~
{: .output}

Something here about what each column means.

We can then see the sequence of these genes by looking at the `.faa` file.
~~~
head bin.6.faa
~~~
{: .bash}

~~~
>JODGFBLK_00001 Lysine--tRNA ligase
MSQEEHNHEELNDQLQVRRDKMNQLRDNGIDPFGARFERTHQSQEVISAYQDLTKEELEE
KAIEVTIAGRMMTKRGKGKAGFAHLQDLEGQIQIYVRKDSVGDDQYEIFKSSDLGDLIGV
TGKVFKTNVGELSVKATSFELLTKALRPLPDKYHGLKDVEQRYRQRYLDLIVNPDSKHTF
IARSKIIQAMRRYLDDHGYLEVETPTMHSIPGGASARPFITHHNALDIPLYMRIAIELHL
KRLIVGGLEKYMNNT
>JODGFBLK_00002 tRNA-dihydrouridine synthase B
MFKIGDIQLKNRVVLAPMAGVCNSAFRLTVKEFGAGLVCAEMVSDKAILYNNARTMGMLY
IDEREKPLSLQIFGGKKETLVEAAKFVDQNTTADIIDINMGCPVPKITKCDAGAKWLLDP
DKIYEMVSAVVDAVDKPVTVKMRMGWDEDHIFAVENAKAVERAGGKAVALHGRTRVQMYE
GTANWDIIKDVKQSVSIPVIGNGDVKTPQDAKRMLDETGVDGVMIGRAALGNPWMIYRTV
QYLETGELKEEPQVREKMAVCKLHLDRLINLKGENVAVREMRKHAAWYLKGVRGNANVRN
EINHCETREEFVQLLDAFTVEVEAKELQNAKVG
>JODGFBLK_00003 hypothetical protein
MEAEIWGRRIRAFRKLKGYTQEGFAKALGISVSILGEIERGNRLPSAAIIQGAADVLNIS
ADELAPPEKDNE
~~~
{: .output}


## Relating these genes to an online database

Introduction to KEGG

* Download the *.faa file.

Once you have dowloaded the Prokka `*.faa` file onto your local computer you can upload this onto [BlastKOALA](https://www.kegg.jp/blastkoala/) in order to annotate the genes with K numbers to relate back to the KEGG database.

<img src="{{ page.root }}/fig/04_03_blastkoala.png" alt="a screenshot of the blastKoala upload page" />

<img src="{{ page.root }}/fig/04_03_blastkoala_upload.png" alt="a screenshot of the blastKoala upload button" />

You should click on the "Choose file" button and navigate to where your `*.faa` file is located on your computer.

We will leave all the options on default, but you need to add in your email address so BlastKOALA can email you to start the job.

Once you have pressed submit you should be re-directed to a screen that says "Request accepted" you will also be sent an email with two links, one to submit the job and one to cancel. **Make sure you press the submit link as your job will not be running without it!** (& make sure to check your spam!)
Once you have pressed the submit link in the email you should be redirected to a BlastKOALA page that says "Job submitted"

As this is an online server your job will sit in a queue with other peoples jobs so it may take a bit of time (hours!) before it is run. You will recieve an email when the job has completed.  

Once the job has been completed you will receive a link by email.  

<img src="{{ page.root }}/fig/04_03_blastkoala_out.png" alt="a screenshot of the blastKoala output" />

From this you can explore the annotated MAG. You can view or download the data and use the KEGG mapper reconstruct pathway option to see how these genes interact with each other. 

## Building a tree from the 16S sequence

Prokka is able to identify any 16S sequences present in our MAGs. This is useful to build a quick taxonomic tree to see what other organisms our MAG relates to.  
We will first search for the presence of 16S seqences in the prokka output.
While still logged into the instance, navigate to the prokka output directory you generated earlier.  
Once in that directory we can use grep to search for 16S in the tsv file with the following command:
~~~
grep 16S *.tsv
~~~
{: .bash}
Depending on the MAG you have run this on you should see similar to the below:
~~~
JODGFBLK_00079	rRNA	1547				16S ribosomal RNA
JODGFBLK_00089	rRNA	1553				16S ribosomal RNA
JODGFBLK_00259	rRNA	1547				16S ribosomal RNA
JODGFBLK_00414	rRNA	583				16S ribosomal RNA (partial)
JODGFBLK_00417	rRNA	1546				16S ribosomal RNA
~~~
{: .output}

Note: if you don't get an output here it may be that your MAG doesn't have any 16S sequences present, which means you may have run Prokka on a less complete MAG. You should double check your output from CheckM and pick a MAG that is highly complete to run through Prokka instead.  

From our output we can see that there are 4 full size 16S ribosomal RNA genes present in our data and one partial one.   

The next step is to pull out the sequence of this 16S rRNA gene and run it through a BLAST database.  
Using the program `seqkit` with the `grep` option (see [seqkit grep](https://bioinf.shenwei.me/seqkit/usage/#grep)) we can pull out the sequence from the `*.ffn` file. (Note we are using the `*.ffn` file here as this will give our 16S sequences in nucleotide format).  

The basic format of this command is the following, `-p` indicates the pattern to search which in our case is the prokka ID. This is the alphanumeric string that is in the first column from the grep output above.
~~~
seqkit grep -p <prokka_id> <prokka.ffn>
~~~
{: .bash}
So in our case this command would be:

~~~
$ seqkit grep -p JODGFBLK_00079 bin.6.ffn
~~~
{: .bash}

You can see the output of our command below:

> ## 16S sequence
> ~~~
> >JODGFBLK_00079 16S ribosomal RNA
> TCGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTC
> GAGCGGACAGATGGGAGCTTGCTCCCTGATGTTAGCGGCGGACGGGTGAGTAACACGTGG
> GTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGCTTGT
> TTGAACCGCATGGTTCAAACATAAAAGGTGGCTTCGGCTACCACTTACAGATGGACCCGC
> GGCGCATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCAACGATGCGTAGCCGACCTG
> AGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAG
> TAGGGAATATTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGG
> TTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTACCGTTCGAATAGGGCGGTA
> CCTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATAC
> GTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTCCTTAA
> GTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAG
> TGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAA
> CACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGG
> AGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAG
> GGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGT
> CGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTT
> TAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCCTAGAGA
> TAGGACGTCCCCTTCGGGGGCAGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTC
> GTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGATCTTAGTTGCCAGCATTC
> AGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAA
> ATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACAGAACAAAGGGCA
> GCGAAACCGCGAGGTTAAGCCAATCCCACAAATCTGTTCTCAGTTCGGATCGCAGTCTGC
> AACTCGACTGCGTGAAGCTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATA
> CGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTC
> GGTGAGGTAACCTTTTTAGGAGCCAGCCGCCGAAGGTGGGACAGATGATTGGGGTGAAGT
> CGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTT
> ~~~
> {: output}
{: .solution}

Now we have the 16S rRNA sequence we can upload this to BLAST and search the 16S database to see what organisms this MAG relates to.

## BLAST

We will be using BLAST (Basic Local Alignment Search Tool) which is an algorithm to find regions of similarity between biological sequences. BLAST is a very popular program in bioinformatics so you may be familiar with the online BLAST server NCBI runs.

We will be using the online server, available at [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

<img src="{{ page.root }}/fig/04_03_blast.png" alt="a screenshot of the blast website" />

Once you have navigated to the page you should select the button that says "Nucleotide BLAST".

<img src="{{ page.root }}/fig/04_03_blast2.png" alt="a screenshot of the submission page" />

Once on this [website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) you should select the rRNA/ITS database under the "Choose Search Set" and then you can paste in the 16S sequence in the box titled "Enter accession number(s), gi(s), or FASTA sequence(s)"

Your screen should look something like the below:

<img src="{{ page.root }}/fig/04_03_blast3.png" alt="a screenshot of the query sequence and 16S/ITS rRNA selected" />

Once you have added in your sequence and selected the rRNA/ITS database you can click the blue BLAST button.
This will then send your job off to the queue, how long it takes from here will depend how busy the NCBI server is - usually it will only take a couple of minutes as the sequence length and the size of the database we're using is small. You should make sure you leave the tab open as this will be where your results arrive.

Once your job has finished you should see an output like below.

<img src="{{ page.root }}/fig/04_03_blast4.png" alt="Output of a BLAST search" />

From here you can explore what sequences have been aligned to your 16S sequence, in the "Descriptions", "Graphic Summary", "Alignments" and "Taxonomy" tabs. You can also browse the "Distance tree of results" to see where your 16S sequence lies in relation to the others.




{% include links.md %}
