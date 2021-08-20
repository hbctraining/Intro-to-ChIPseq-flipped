---
title: "Visualization of peaks"
author: "Meeta Mistry, Jihe Liu"
date: "Aug 19th, 2021"
---

Approximate time:

**Link to issue describing the modifications to be made:** https://github.com/hbctraining/Intro-to-ChIPseq-flipped/issues/11

## Learning Objectives
* Understand what is bigWig file
* Learn how to generate bigWig files using `deepTools`

## Visualization of ChIP-seq data

The first stage of ChIP-seq analysis uses common processing pipelines, which involves the alignment of raw reads to the genome, data filtering, and identification of enriched signal regions (peak calling). In the second stage, individual software programs allow detailed analysis of those peaks, biological interpretation, and visualization of ChIP-seq results.

There are various strategies for visualizing enrichment patterns, and we will explore a few of them. To start, we will create bigWig files for our samples, a standard file format commonly used for ChIP-seq data visualization.

## Creating bigWig files

The general procedure is to take our alignment files (BAM) and convert them into bigWig files. The bigWig format is an indexed binary format useful for dense, continuous data to be displayed in a genome browser as a graph/track. The bigWig files could also be used as inputs for some visualization commands we will be running in `deepTools`. 

[`deepTools`](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html) is a suite of Python tools developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq, or MNase-seq. `deepTools` has a wide variety of commands that go beyond what we will cover in this lesson. We encourage you to look through the docuementation and explore on your own time. 

<p align="center">
<img src="../img/bam_to_bigwig.png" width="700">
</p>

*Image acquired from the [deepTools documentation](http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html?highlight=bigwig)*

Start an interactive session with 6 cores. *If you have already logged on to a compute node, you will want to exit and start a new session*.

```bash
$ srun --pty -p interactive -t 0-12:00 --mem 8G -c 6 --reservation=HBC2 bash
```

We will begin by creating a directory for the visualization output and loading the required modules to run `deepTools`.

```bash
$ cd ~/chipseq_workshop/results/
$ mkdir -p visualization/bigWig
```

```bash
$ module load gcc/6.2.0  python/2.7.12
$ module load deeptools/3.0.2 
```

We then need to **create an index file for the BAM file**. To perform some functions on the BAM file, many tools require an index. Think of an index located at the back of a textbook - when you are interested in a particular subject, you look for the keyword in the index and identify the pages that contain the relevant information. Similarily, indexing the BAM file aims to achieve fast retrieval of alignments overlapping a specified region without going through the whole alignment file. 

We will use [SAMtools](http://samtools.sourceforge.net/) again, specifically the **`samtools index`** command, to index a BAM file.

Let's load the `samtools` module:

```bash
$ module load samtools/1.9
```

Create index for the `wt_sample2_chip_final` BAM file that we created in earlier lesson:

```bash
$ samtools index ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam
```

There are two useful commands in `deepTools` to create bigWig files: `bamCoverage` and `bamCompare`. The former takes in a single BAM file and return a bigWig file. The latter normalizes two files to each other (i.e. ChIP sample relative to input), and returns a single bigWig file.

We will use the `bamCoverage` command to **create a bigWig file for `wt_sample2_chip`**. We will specify `binSize` of 20, as an additional parameter. There are a few other parameters that you could explore (but we will not use). 

* `normalizeUsing`: Possible choices: RPKM, CPM, BPM, RPGC.
* `binSize`: size of bins in bases
* `smoothLength`: defines a window, larger than the `binSize`, to average the number of reads over. This helps produce a more continuous plot.
* `centerReads`: reads are centered with respect to the fragment length as specified by `extendReads`. This option is useful to get a sharper signal around enriched regions.

```bash
$ bamCoverage -b ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20
```

Alternatively, we could use `bamCompare` to **create a bigWig file in which we normalize the ChIP against the input**. The command is quite similar to `bamCoverage`, except that it require two files as input (`b1` and `b2`). [Here](https://deeptools.readthedocs.io/en/develop/content/help_faq.html#when-should-i-use-bamcoverage-or-bamcompare) are details about the difference between `bamCompare` and `bamCoverage`.

```bash
## DO NOT RUN

$ bamCompare -b1 ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-b2 ~/chipseq_workshop/results/bowtie2/wt_sample2_input_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20
```

