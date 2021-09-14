---
title: "Handling peak files with bedtools"
author: "Meeta Mistry, Jihe Liu"
date: "Aug 17th, 2021"
---

Contributors: Meeta Mistry, Radhika Khetani, Jihe Liu

Approximate time:

## Learning Objectives

* Understand the BED file format and related file formats
* Utilize the bedtools suite of tools filter peak calls 
* Utilize the bedtools suite to compare peaks betweeen replicates


## Handling peak calls 

In this lesson, we will introduce you to an important file format that you will encounter when working with peak calls called the **BED format**. We will discuss the output files that we obtained from MACS2 peak calling, specifically describing the contents of the narrowPeak files and how it relates to BED. You will then get acquainted with **`bedtools`, a new suite of tools that is very helpful when working with BED files and other related file formats**, and use it to complete the folowing tasks with the WT and KO peak calls from this PRDM16 dataset:

1. Filter out peaks that overlap with the blacklisted regions
2. Assess the replicate concordance within sample groups, to see how many peaks are reproducible. 


### BED file formats

#### BED
Take content from slide deck here

#### narrowPeak
A narrowPeak (.narrowPeak) file is used by the ENCODE project to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED 6+4 format, which means the first 6 columns of a standard BED file  with **4 additional fields**:

**DON'T HAVE TO KEEP THIS IMAGE** can write out the text version or whatever works best 

<p align="center">
<img src="../img/narrowPeak.png">
</p>

### `bedtools`

The **bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks**. The general idea is that genome coordinate information can be used to perform relatively simple arithmetic, like combining, subsetting, intersecting etc., to obtain desired information. [bedtools](http://bedtools.readthedocs.org/en/latest/index.html) was devloped in [Aaron Quinlan's group](http://quinlanlab.org/) at University of Utah, and is widely used amongst the bioinformatics community. It is such an easy and versatile tool to perform these tasks described above. 

<p align="center">
<img src="../img/bedtools.png" width="700">
</p>

As the name implies, this suite of tools works with **Bed** files, but it also works with other file formats that have genome coordinate information. 

<p align="center">
<img src="../img/bedtools-basic.png" width="600">
</p>

> **NOTE:** When working with multiple files to perform arithmetic on genomic coordinates, it is essential that all files have coordinate information from the same version of the genome and the coordinate system (0-based or 1-based)!

`bedtools` is available as a module on O2. To set yourself up for the rest of the lesson, make sure you are **logged into O2 and on a compute node**.

```bash
$ srun --pty -p interactive -t 0-12:00 --mem 2G /bin/bash
```

Next, load the modules for `bedtools` and `samtools`:
	
```bash
$ module load gcc/6.2.0 bedtools/2.26.0 samtools/1.3.1
```


### `bedtools intersect`

The [`bedtools intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) command within bedtools evaluates A (file 1) and finds regions that overlap in B (file 2). We will use this command to do both the filtering of peaks (from blacklisted regions) and assessing the overlap of peaks (between replicates).

<p align="center">
<img src="../img/bedtools_intersect.png" width="600">
</p>


To find out more information on the parameters available when intersecting, use the help flag:

```bash
$ bedtools intersect -h
```

Alternatively, you can use the [web-based documentation](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) which is much easier to read. The "options" summary shows the **large number of arguments available that allow us to do so much more than just a simple A versus B comparison**. It allows one to have fine control as to how the intersections are reported, and work with different types of files - amongst may other features.

### Filtering peaks overlapping with blacklist regions

We discussed blacklisted regions in the [filtering lesson](05_filtering_BAM_files.md), as it is commonplace to filter before peak calling. When it is performed at the stage the `bedtools intersect` is also used, the difference being the input file type (BAM instead of BED). The filtering works just as well if applied post-peak calling.

The blacklisted regions typically appear uniquely mappable so simple mappability filters do not remove them. These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats.

<p align="center">
<img src="../img/blacklist.png" width="600">
</p>

We have a BED file of blacklist regions for mouse `mm10` prepared at `/n/groups/hbctraining/harwell-datasets/workshop_material/reference/mm10-blacklist.v2.bed`. Copy this file over to your project into the `reference_data` folder.

```bash
cp /n/groups/hbctraining/harwell-datasets/workshop_material/reference/mm10-blacklist.v2.bed  ~/chipseq_workshop/reference_data
```

> **How were the 'blacklists compiled?** These blacklists were empirically derived from large compendia of data using a combination of automated heuristics and manual curation. Blacklists were generated for various species and genome versions including human, mouse, worm and fly. The lists can be ]downloaded here](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/) or from [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). For human, they used 80 open chromatin tracks (DNase and FAIRE datasets) and 12 ChIP-seq input/control tracks spanning ~60 cell lines in total. These blacklists are applicable to functional genomic data based on short-read sequencing (20-100bp reads). These are not directly applicable to RNA-seq or any other transcriptome data types.

Next, we will want to navigate to our results directory:

```bash
$ cd ~/chipseq_workshop/results/
```


To filter out blacklisted region for `wt_sample1`, we use the following code. Note that we use the flag `-v`, to report entries in A that have no overlaps with B.


More information about the blacklist region is described in this [paper](https://www.nature.com/articles/s41598-019-45839-z), and we downloaded the blacklist regions from [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). To filter out blacklisted region for `wt_sample1`, we use the following code. Note that we use the flag `-v`, to report entries in A that have no overlaps with B.

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



