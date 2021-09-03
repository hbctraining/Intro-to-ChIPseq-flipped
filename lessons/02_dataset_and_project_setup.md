---
title: "Introduction to the dataset and project organization"
author: "Mary Piper, Radhika Khetani, Meeta Mistry"
date: "Sept 1st, 2021"
---

Approximate time: 45 minutes

## Learning Objectives

* Introducing the dataset and biological context
* Recognize the need for data management and project organization


## Introduction to the dataset

For this workshop we will be working with ChIP-seq data from a recent publication in Neuron by Baizabal et al. (2018) [[1]](https://doi.org/10.1016/j.neuron.2018.04.033). 

The authors in this paper sought to understand how chromatin-modifying enzymes function in neural stem cells to establish the epigenetic landscape that determines cell type and stage-specific gene expression. Chromatin-modifying enzymes are transcriptional regulators that control gene expression through covalent modification of DNA or histones. 

<p align="center">
<img src="../img/05Epigenetics06-edit.png" width="400">
</p>

_Image adapted from: [American Society of Hematology](https://www.hematology.org/research/ash-agenda-for-hematology-research/epigenetic-mechanisms)_

### PRDM16

The transcriptional regulator **PRDM16 is a chromatin-modifying enzyme** that belongs to the larger PRDM (Positive Regulatory Domain) protein family, that is structurally defined by the **presence of a conserved N-terminal histone methyltransferase PR domain** (Hohenauer and Moore, 2012). 

* PRDM16 has been shown to function in vitro as a histone 3 lysine 9 (H3K9) and histone 3 lysine 4 (H3K4) mono-methyltransferase (Pinheiro et al., 2012, Zhou et al., 2016). 
* PRDM16 also regulates gene expression by forming complexes with transcriptional co-factors and other histone-modifying proteins (Chi and Cohen, 2016). 
* PRDM16 was previously shown to control embryonic and post-natal neural stem cell maintenance and differentiation in the brain (Chuikov et al., 2010, Inoue et al., 2017, Shimada et al., 2017). 


**How PRDM16 functions to regulate transcriptional programs in the developing cerebral cortex remains largely unknown.**

In this paper, the authors use various techniques to identify and validate the targets and activities of PRDM16, including ChIP-seq, bulk RNA-seq, FACS, in-situ hybridization and immunofluorescent microscopy on brain samples from embryonic mice and a generation of PRDM16 conditional knockout mice. 

<p align="center">
<img src="../img/graphical_abstract.png" width="500">
</p>


From the RNA-seq data, they found that the absence of PRDM16 in cortical neurons resulted in the misregulation of over a thousand genes during neurogenesis. To identify the subset of genes that are transcriptional targets of PRDM16 and to understand how these genes are directly regulated, they performed chromatin immunoprecipitation followed by sequencing (ChIP-seq).


**Hypothesis:** Does the histone methyltransferase PRDM16 works with enhancer elements to either silence or activate expression of sets of genes that impact the organization of the cerebral cortex?


### Raw data

From this study we are using the [ChIP-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111657) data which is publicly available in the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra?term=SRP134733).

> **NOTE:** If you are interested in how to obtain publicly available sequence data from the SRA we have some materials on this [linked here](https://hbctraining.github.io/Accessing_public_genomic_data/lessons/downloading_from_SRA.html).

### Metadata

In addition to the raw sequence data we also need to collect **information about the data**, also known as **metadata**.  We are usually quick to want to begin analysis of the sequence data (FASTQ files), but how useful is it if we know nothing about the samples that this sequence data originated from? 

Some relevant metadata for our dataset is provided below:

* Whole brain lysates were obtained from **mice at E15.5** (when upper layer neurons are being generated)
* Approximately **20-40 million cortical cells** were isolated for each sample
* After the ChIP protocol, genomic DNA was purified, end repaired, ligated with barcoded adaptors, amplified for 11 PCR cycles
* Library **fragments in the range of 100-800 bp** were size-selected using agarose gel electrophoresis followed by DNA gel extraction
* Libraries were sequenced in an **Illumina HiSeq 2500** sequencer to a sequencing depth of **30-40 million reads per sample**.

> All of the above pertains to both WT and Prdm16 conditional knock-out (cKO) mouse (Emx1Ires-Cre; Prdm16flox/flox)

Our dataset consists of two WT samples and two cKO samples. For each of the IP samples, we have a corresponding input sample as illustrated in the schematic below.

**Add image here**

## Implementing data management best practices

In a [previous lesson](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/04a_data_organization.html) we describe the data lifecycle and the **different aspects to consider when working on your own projects**. Here, we implement some of those strategies to get ourselves setup before we begin with any analysis.

<p align="center">
<img src="../img/data-lifecycle-base.png" width="500">
</p>

_Image acquired from the [Harvard Biomedical Data Management Website](https://datamanagement.hms.harvard.edu/data-lifecycle)_

### Planning and organization

For each experiment you work on and analyze data for, it is considered best practice to get organized by creating a planned storage space (directory structure). We will start by creating a directory that we can use for the rest of the workshop. First, make sure that you are in your home directory.

```bash
$ cd
$ pwd
```

This should return `/home/rc_training`. Create the directory `rnaseq` and move into it.

```bash
$ mkdir rnaseq
$ cd rnaseq
```

## Project Organization

Since we are going to be working with this data on our remote server, **O2**, we first need to log onto the server. 

Type in the following command with your username to login:

```bash
ssh username@o2.hms.harvard.edu
```

Next we will start an interactive session on O2 with 2 cores (add the `-n 2`):

```bash
$ srun --pty -p interactive -t 0-12:00 --mem 1G -n 1 --reservation=HBC2 /bin/bash
```

**Make sure that your command prompt is now preceded by a character string that contains the word "compute".**

>**NOTE:** We are using the `--reservation` argument during class since we have a dedicated set of computers reserved so that commands run quickly. When starting an interactive session outside of class you will need to leave out this argument. You may also want to increase othe resources (i.e. memory and number of cores) depending on what you plan on doing.
>
>```bash
>$ srun --pty -p interactive -t 0-12:00 --mem 8G -n 2 /bin/bash
>```


## Data Management

One of the most important parts of research that involves large amounts of data, is how best to manage it. We tend to prioritize the analysis, but there are many other important aspects that are  often overlooked in the excitement to get a first look at new data. 

The data management lifecycle displayed below, courtesy of the [HMS Data Management Working Group](https://datamanagement.hms.harvard.edu/hms-data-management-working-group), illustrates some things to consider beyond the data creation and analysis components:

<img src="../img/data_life_cycle_gouldv2.png" width="350">

_Image aquired from the [Harvard Biomedical Data Management Website](https://datamanagement.hms.harvard.edu/hms-data-lifecycle)_

We will cover some parts of this lifecycle by talking about best practices for the **Research** half of the above lifecycle. Later in this workshop we will talk a little more about the data storage. For more information about the full lifecycle and more guidelines for data management, please look at the resources linked below.

**Resources**

* The [HMS Data Management Working Group's website](https://datamanagement.hms.harvard.edu/)
* A guide from the [Harvard library](http://guides.library.harvard.edu/dmp).
* Sign-up for the [DMWG quarterly newsletter](https://harvard.us13.list-manage.com/subscribe?u=d3fee19ad91470512cdd564dd&id=13642f2d02) for helpful tips, classes and events related to data management


### Planning

You should approach your sequencing project in a very similar way to how you do a biological experiment, and ideally, begins with **experimental design**. We're going to assume that you've already designed a beautiful sequencing experiment to address your biological question, collected appropriate samples, and that you have enough statistical power.

During this stage it is important to keep track of how the experiment was performed and clearly tracking the source of starting materials and kits used. It is also best practice to include information about any small variations within the experiment or variation relative to standard experiments. 

### Organization

Every computational analysis you do is going to spawn many files, and inevitability you'll want to run some of those analyses again. For each experiment you work on and analyze data for, it is considered best practice to get organized by creating a planned storage space (directory structure).

We will start by creating a directory that we can use for the rest of the ChIP-seq session.

First, make sure that you are in your home directory.

```bash
$ cd
$ pwd
```
This should return `/home/username`.

Create a `chipseq` directory and change directories into it:

```bash
$ mkdir chipseq

$ cd chipseq
```

Now that we have a project directory, we can set up the following structure within it to keep files organized.

```bash
chipseq/
├── logs/
├── meta/
├── raw_data/
├── reference_data/
├── results/
│   ├── bowtie2/
│   └── fastqc/
└── scripts/
```

```bash
$ mkdir raw_data reference_data scripts logs meta

$ mkdir -p results/fastqc results/bowtie2

$ tree     # this will show you the directory structure you just created
```

> **NOTE:** We are using the parents flag (`-p` or `--parents`) with `mkdir` to complete the file path by creating any parent directories that do not exist. In our case, we have not yet created the `results` directory and so since it does not exist it will be created. This flag can be very useful when scripting workflows. 

**This is a generic directory structure and can be tweaked based on personal preference and analysis workflow.**

- `logs`: to keep track of the commands run and the specific parameters used, but also to have a record of any standard output that is generated while running the command. 
- `meta`: for any information that describes the samples you are using, which we refer to as [metadata](https://datamanagement.hms.harvard.edu/metadata-overview). We will discuss this in more detail as it pertains to our example dataset, later in this lesson.
- `raw_data`: for any **unmodified** (raw) data obtained prior to computational analysis here, e.g. FASTQ files from the sequencing center. We strongly recommend leaving this directory unmodified through the analysis.
- `reference_data`: for known information related to the reference genome that will be used in the analysis, e.g. genome sequence (FASTA), gene annotation file (GTF) associated with the genome.
- `results`: for output from the different tools you implement in your workflow. Create sub-folders specific to each tool/step of the workflow within this folder. 
- `scripts`: for scripts that you write and use to run analyses/workflow.


Now that we have the directory structure created, let's copy over the data to perform our quality control and alignment, including our FASTQ files and reference data files:

```bash
$ cp /n/groups/hbctraining/chip-seq/raw_fastq/*fastq raw_data/

$ cp /n/groups/hbctraining/chip-seq/reference_data/chr12* reference_data/
```

Now we are all set up for our analysis!

> #### File naming conventions
> 
> Another aspect of staying organized is making sure that all the filenames in an analysis are as consistent as possible, and are not things like `alignment1.bam`, but more like `20170823_kd_rep1_gmap-1.4.bam`. [This link](https://datamanagement.hms.harvard.edu/file-naming-conventions) and [this slideshow](http://www2.stat.duke.edu/~rcs46/lectures_2015/01-markdown-git/slides/naming-slides/naming-slides.pdf) have some good guidelines for file naming dos and don'ts.


### Documentation

**Documentation doesn't stop at the sequencer!** Keeping notes on what happened in what order, and what was done, is essential for reproducible research.

#### Log files

In your lab notebook, you likely keep track of the different reagents and kits used for a specific protocol. Similarly, recording information about the tools and parameters is important for documenting your computational experiments. 

- **Make note of the software you use.** Do your research and find out what tools are best for the data you are working with. Don't just work with tools that you are able to easily install.
- **Keep track of software versions.** Keep up with the literature and make sure you are using the most up-to-date versions.
- **Record information on parameters used and summary statistics** at every step (e.g., how many adapters were removed, how many reads did not align)
    - A general rule of thumb is to test on a single sample or a subset of the data before running your entire dataset through. This will allow you to debug quicker and give you a chance to also get a feel for the tool and the different parameters.
    - Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture. You can redirect standard output with the `>` symbol which is equivalent to `1> (standard out)`; other tools might require you to use `2>` to re-direct the `standard error` instead.
    
#### README files

After setting up the directory structure and when the analysis is running it is useful to have a **[README file](https://datamanagement.hms.harvard.edu/readme-files) within your project directory**. This file will usually contain a quick one line summary about the project and any other lines that follow will describe the files/directories found within it. An example README is shown below. Within each sub-directory you can also include README files to describe the analysis and the files that were generated.

```
## README ##
## This directory contains data generated during the Intro to ChIP-seq course
## Date: 

There are six subdirectories in this directory:

raw_data : contains raw data
meta:  contains...
logs:
reference_data:
results:
scripts:
```

*** 

### Homework Exercise

- Create a README for the `chipseq/` folder (hint: use `vim` to create the file). Give a short description of the project and as homework add brief descriptions of the types of files you will be storing within each of the sub-directories. 

***




***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
