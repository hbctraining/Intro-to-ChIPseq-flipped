---
title: "Handling replicates"
author: "Meeta Mistry, Jihe Liu"
date: "Aug 17th, 2021"
---

Contributors: Meeta Mistry, Radhika Khetani, Jihe Liu

Approximate time:

**Link to issue describing the modifications to be made:** https://github.com/hbctraining/Intro-to-ChIPseq-flipped/issues/10

## Learning Objectives

* Learn how to filter out peaks at blacklisted regions
* Learn how to find overlap peaks betweeen replicates


## Filtering peaks

In this section, our goal is to refine our called peaks from macs2. We will do two consecutive filterings: first, filter out peaks that are within the blacklisted region; second, find overlap peaks between two biological replicates. We will use a suite of tools called `bedtools` to perform the task.

### `bedtools`

The general idea is that genome coordinate information can be used to perform relatively simple arithmetic, like combining, subsetting, intersecting etc., to obtain desired information. [bedtools](http://bedtools.readthedocs.org/en/latest/index.html) from [Aaron Quinlan's group](http://quinlanlab.org/) at University of Utah is such an easy and versatile tool to perform these tasks. 

<p align="center">
<img src="../img/bedtools.png" width="700">
</p>

As the name implies, this suite of tools works with **Bed** files, but it also works with other file formats that have genome coordinate information. 

<p align="center">
<img src="../img/bedtools-basic.png" width="600">
</p>

> **NOTE:** When working with multiple files to perform arithmetic on genomic coordinates, it is essential that all files have coordinate information from the same version of the genome and the coordinate system (0-based or 1-based)!

### Setting up

Let's start an interactive session. 

```bash
$ srun --pty -p short -t 0-12:00 --mem 8G --reservation=HBC bash	

$ cd ~/chipseq_workshop/results/
```
	
Load the modules for `bedtools` and `samtools`:
	
```bash
$ module load gcc/6.2.0 bedtools/2.26.0 samtools/1.3.1
```

### `bedtools intersect`

The [`bedtools intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) command within bedtools evaluates A (file 1) and finds regions that overlap in B (file 2). We will use this command to do the filtering.

<p align="center">
<img src="../img/bedtools_intersect.png" width="600">
</p>

To find out more information on the parameters available when intersecting, use the help flag:

```bash
$ bedtools intersect -h
```

### filtering out peaks in blacklisted regions

We discussed blacklisted regions in the previous lesson, where we mentioned that although we could have filtered out blacklisted regions from BAM files, we would instead perform the filtering after peak calling. So here we are! As you will see, filtering blacklisted regions using `bedtools` is quick and straightforward.

We have deposited the blacklist regions for mouse `mm10` version in `~/chipseq_workshop/reference/mm10-blacklist.v2.bed`. More information about the blacklist region is described in this [paper](https://www.nature.com/articles/s41598-019-45839-z), and we downloaded the blacklist regions from [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). To filter out blacklisted region for `wt_sample1`, we use the following code. Note that we use the flag `-v`, to report entries in A that have no overlaps with B.

```bash
bedtools intersect \
-v \
-a ~/chipseq_workshop/results/macs2/wt_sample1_peaks.narrowPeak \
-b ~/chipseq_workshop/reference/mm10-blacklist.v2.bed \
> ~/chipseq_workshop/results/macs2/wt_sample1_peaks_filtered.bed
```

Similarly, we could filter out blacklisted regions for `wt_sample2`:

```bash
bedtools intersect \
-v \
-a ~/chipseq_workshop/results/macs2/wt_sample2_peaks.narrowPeak \
-b ~/chipseq_workshop/reference/mm10-blacklist.v2.bed \
> ~/chipseq_workshop/results/macs2/wt_sample2_peaks_filtered.bed
```

> **NOTE:** The narrowPeak file also follows the bed file format. That's why we could use it as an input, even though it does not end with the suffix `.bed`.

We could use `wc -l` command to check how many peaks are filtered out because they are located at the blacklisted region:

```bash
# Number of peaks before the filtering
wc -l ~/chipseq_workshop/results/macs2/wt_sample2_peaks.narrowPeak

# Number of peaks after the filtering
wc -l ~/chipseq_workshop/results/macs2/wt_sample2_peaks_filtered.bed
```

### Finding overlapping peaks between replicates

Next, as we have two replicates per condition, we would like to find out the overlapping peaks between replicates - these will be confident peaks we identify as the final result. Similar to the previous step, we could use the `bedtools intersect` command again. The only difference is that we will use three different flags here:

- `-wo`: Write the original A (file 1) and B (file 2) entries plus the number of base pairs of overlap between the two features.
- `-f`: Minimum overlap required as a fraction of A. The value ranges from 0 to 1. We will use 0.3, requiring the overlap region being at least 30% of A.
- `-r`: Require that the fraction overlap be reciprocal for A and B. Together with the `-f` flag above, we require the overlap region being at least 30% of B as well.

The below code generates the overlapping peaks: 

```bash
$ bedtools intersect \
-wo -f 0.3 -r \
-a ~/chipseq_workshop/results/macs2/wt_sample1_peaks_filtered.bed \
-b ~/chipseq_workshop/results/macs2/wt_sample2_peaks_filtered.bed \
> ~/chipseq_workshop/results/macs2/wt_peaks_final.bed
```

Finally, let's check how many confident peaks are present now:

```bash
wc -l ~/chipseq_workshop/results/macs2/wt_peaks_final.bed
```

> **_Historical Note_:** A simpler heuristic for establishing reproducibility was previously used as a standard for depositing ENCODE data and was in effect when much of the currently available data was submitted. According to this standard, either 80% of the top 40% of the peaks identified from one replicate using an acceptable scoring method should overlap the list of peaks from the other replicate, OR peak lists scored using all available reads from each replicate should share more than 75% of regions in common. As with the current standards, this was developed based on experience with accumulated ENCODE ChIP-seq data, albeit with a much smaller sample size.


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*



