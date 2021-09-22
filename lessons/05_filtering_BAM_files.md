---
title: "Alignment and filtering"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu"
date: "Aug 10th, 2021"
---

Contributors: Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu

Approximate time:

## Learning Objectives

* Understand the purpose of filtering alignment reads
* Learn to perform filtering with sambamba and samtools

## Filtering reads

<p align="center">
 <img src="../img/chipseq_filterworkflow_sept2021.png" width="600">
</p>


A key issue when working with a ChIP-seq data is to **move forward with only the uniquely mapping reads**.  Allowing for multiply mapped reads increases the number of usable reads and the sensitivity of peak detection; however, the number of false positives may also increase [[1]](https://www.ncbi.nlm.nih.gov/pubmed/21779159/). To increase our confidence in peak calling and improve data reproducibility, we need to **filter out both multi-mapping reads and duplicate reads**.

* Multi-mapping reads are reads that are mapping to multiple loci on the reference genome.
* Duplicate reads are reads that map at the exact same location, with the same coordination and the same strand. These duplicates can arise from experimental artifacts, but can also contribute to genuine ChIP-signal.
    * **The bad kind of duplicates:** If initial starting material is low, this can lead to overamplification of this material before sequencing. Any biases in PCR will compound this problem and can lead to artificially enriched regions. 
    * **The good kind of duplicates:** You can expect some biological duplicates with ChIP-seq since you are only sequencing a small part of the genome. This number can increase if your depth of coverage is excessive or if your protein only binds to few sites. If there are a good proportion of biological dupicates, removal can lead to an underestimation of the ChIP signal. 

> #### Some additional notes on duplicates
> Most peak calling algorithms also implement methods to deal with duplicate reads. While they are commonly removed prior to peak calling, another option is to leave them now and deal with them later. **Skip the duplicate filtering at this step if**:
> * You are planning on performing a differential binding analysis.
> * You are expecting binding in repetitive regions (also, use paired-end sequencing) 
> * You have included UMIs into your experimental setup.

The older version of Bowtie2 had an argument that allowed us to easily perform filtering during the alignment process. We do not have this option with Bowtie2 and so the filtering will be done with the use of a tool called [sambamba](https://lomereiter.github.io/sambamba/). Sambamba is an open source tool that provides methods for working with SAM/BAM files, similar to samtools, except with faster processing times and in some cases added functionality. 

This **lesson will consist of two steps**:

1. Sort BAM files by genomic coordinates (using `samtools`).
2. Filter the reads to keep only uniquely mapping reads (using `sambamba`). This will also remove any unmapped reads.

Before we begin, you will want to make sure you are **logged into O2.** To start an interactive session with 2 cores and 10G of memory (sorting can be memory-intensive) us the command below:

```bash
$ srun --pty -p interactive -t 0-2:30 --mem 10G -n 1 --reservation=HBC2 /bin/bash
```
> **Make sure that your command prompt is now preceded by a character string that contains the word "compute".**

We will also load the required modules for this lesson:

```bash
module load gcc/6.2.0 samtools/1.13 sambamba/0.7.1
```


### 1. Sort BAM files by genomic coordinates

Before we can do the filtering, we need to sort our BAM alignment files by genomic coordinates (instead of by name). To perform the sorting, we will use [Samtools](http://www.htslib.org/), a tool we previously used when coverting our SAM file to a BAM file. 

The command we use this time is `samtools sort` with the following parameter:

* `-o`: /path/to/output/file

Move into the `bowtie2` directory and then run the command. Once complete, you should see a new file generated with the output file name your provided.

```bash
$ cd ~/chipseq_workshop/results/bowtie2/
$ samtools sort wt_sample2_chip.bam -o wt_sample2_chip_sorted.bam
```

> **NOTE**: You will need the BAM file generated from the [alignment lesson](04_alignment_using_bowtie2.md). If you do not have this file, please copy over the BAM file to your directory:
>
> ```bash
> $ cp /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip.bam ~/chipseq_workshop/results/bowtie2/wt_sample2_chip.bam
> ``` 


### 2. Filter the reads to keep only uniquely mapping reads

Next, we can filter the sorted BAM files to keep only uniquely mapping reads. We use the `sambamba view` command with the following parameters:

* `-t`: number of threads(cores)
* `-h`: print SAM header before reads
* `-f`: format of output file (default is SAM)
* `-F`: set [custom filter](https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax) - we will be using the filter to remove duplicates, multimappers and unmapped reads.

```bash
$ sambamba view -h -t 2 -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
wt_sample2_chip_sorted.bam > wt_sample2_chip_final.bam
```

We filter out unmapped reads by specifying in the filter `not unmapped`, and duplicates with `not duplicate`. Also, among the reads that are aligned, we filter out multimappers by specifying `[XS] == null`. 'XS' is a tag generated by Bowtie2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned and more than one alignment is found for the read.

Now that the alignment files contain only uniquely mapping reads, we are ready to perform peak calling!


> ### Filtering out Blacklisted Regions
> Although we do not perform this step, it is common practice to apply an additional level of filtering to our BAM files. That is, we remove alignments that occur with defined Blacklisted Regions. **We will filter out blacklist regions post-peak calling.**
> 
> Blacklisted regions represent artifact regions that tend to show artificially high signal (excessive unstructured anomalous reads mapping). These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats and typically appear uniquely mappable so simple mappability filters applied above do not remove them. The ENCODE and modENCODE consortia have compiled blacklists for various species and genome versions including human, mouse, worm and fly. These blacklisted regions (coordinate files) can be filtered out from our alignment files before proceeding to peak calling.
> 
> If we wanted to filter blacklist regions at this point in our workflow, we would use the following code:
> 
> ``` 
> # DO NOT RUN
> $ bedtools intersect -v -abam wt_sample2_chip_final.bam -b mm10-blacklist.v2.bed > wt_sample2_chip_final_blacklist_filtered.bam
> ```
> 
> _bedtools is a suite of tools that we will discuss in more detail in a later lesson when blacklist filtering is applied._


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

