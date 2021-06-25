---
title: "Visualization of peaks"
author: "Meeta Mistry"
date: "Thursday July 29th, 2017"
---

Approximate time: 80 minutes

## Learning Objectives

* Visualizing enrichment patterns at particular locations in the genome

## Visualization of ChIP-seq data

The first part of ChIP-sequencing analysis uses common processing pipelines, which involves the alignment of raw reads to the genome, data filtering, and identification of enriched signal regions (peak calling). In the second stage, individual software programs allow detailed analysis of those peaks, biological interpretation, and visualization of ChIP-seq results.

There are various strategies for visualizing enrichment patterns and we will explore a few of them. To start, we will create bigWig files for our samples, a standard file format commonly used for ChIP-seq data visualization.


## Profile plots and heatmaps

Because many cis-regulatory elements are close to TSSs of their targets, a common visualization technique is to use bigWig files to obtain a global evaluation of enrichment around the TSS. In our example, we will assess enrichment around the TSS and plot this separately for the Nanog and Pou5f1 samples (two replicates in each plot). 

Rather than looking at the TSS for all known genes, we will only look be looking at genes on chromosome 12 in the interest of time. Copy over the BED file which contains the coordinates for all genes on chromosome 12 to the visualization folder.

```bash
$ cp /n/groups/hbctraining/chip-seq/deepTools/chr12_genes.bed ~/chipseq/results/visualization/
```

Before we start plotting our data, we first need to prepare an intermediate file that can be used with the `plotHeatmap` and `plotProfile` commands.

<img src="../img/computeMatrix_overview.png" width="700">


The `computeMatrix` command accepts multiple bigWig files and multiple region files (BED format) to create a count matrix which is the intermediate file. It can also be used to filter and sort regions according to their score. Our region file will be the BED file we just copied over and our bigWig files will be those generated from the full dataset that we have provided for you. Additionally, we will specify a window of +/- 1000bp around the TSS of genes (`-b` and `-a`). For each window, `computeMatrix` will calculate scores based on the read density values in the bigWig files.

First, let's create a matrix for one of the Nanog replicates:

```bash

$ computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R ~/chipseq/results/visualization/chr12_genes.bed \
-S /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Nanog*.bw \
--skipZeros \
-o ~/chipseq/results/visualization/matrixNanog_TSS_chr12.gz \
-p 6 \
--outFileSortedRegions ~/chipseq/results/visualization/regions_TSS_chr12.bed

```

> **NOTE:** Typically, the genome regions are genes, and can be obtained from the [UCSC table browser](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables). Alternatively, you could look at other regions of interest that are not genomic feature related (i.e. binding regions from another protein of interest).

Now, let's create another matrix for the Pou5f1 replicates:

```bash

$ computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R ~/chipseq/results/visualization/chr12_genes.bed \
-S /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Pou5f1*.bw \
--skipZeros \
-p 6 \
-o ~/chipseq/results/visualization/matrixPou5f1_TSS_chr12.gz \
--outFileSortedRegions ~/chipseq/results/visualization/regionsPou5f1_TSS_chr12.bed

```

Using that matrix we can create a **profile plot** which is essentially a density plot that evaluates read density across all transcription start sites. For Nanog, we can see that **Replicate 2 has a particularly higher amount of signal at the TSS compared to Replicate 1**. 

```bash
$ plotProfile -m visualization/matrixNanog_TSS_chr12.gz \
-out visualization/figures/TSS_Nanog_profile.png \
--perGroup \
--colors green purple \
--plotTitle "" --samplesLabel "Rep1" "Rep2" \
--refPointLabel "TSS" \
-T "Nanog read density" \
-z ""

```

<img src="../img/TSS_Nanog_profile.png" width="500">

Alternatively, we could use a **heatmap** to evaluate the same matrix of information:

```bash
$ plotHeatmap -m visualization/matrixNanog_TSS_chr12.gz \
-out visualization/figures/TSS_Nanog_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -4 --zMax 4  
```
<img src="../img/TSS_Nanog_heatmap.png" width="400">


Similarly we can do the same for **Pou5f1. Here, we find that Replicate 1 exhibits stronger signal**.

```bash
$ plotProfile -m visualization/matrixPou5f1_TSS_chr12.gz \
-out visualization/figures/TSS_Pou5f1_profile.png \
--perGroup --colors green purple \
--plotTitle "" --samplesLabel "Rep1" "Rep2" \
--refPointLabel "TSS" -T "Pou5f1 read density" -z ""
```

<img src="../img/TSS_Pou5f1_profile.png" width="400">

```bash
$ plotHeatmap -m visualization/matrixPou5f1_TSS_chr12.gz \
-out visualization/figures/TSS_Pou5f1_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -2 --zMax 2  
```

<img src="../img/TSS_Pou5f1_heatmap.png" width="400">

If we wanted **both images in one single plot**, we can do that with `plotHeatmap` and just removing the `--whatToShow` parameter.

```bash
$ plotHeatmap -m visualization/matrixPou5f1_TSS_chr12.gz \
-out visualization/figures/TSS_Pou5f1_profile-heatmap.png \
--colorMap RdBu \
--zMin -2 --zMax 2  
```

<img src="../img/TSS_Pou5f1_heatmap_and_profile.png" width="400">

> **NOTE:** Both `plotProfile` and `plotHeatmap` have many options, including the ability to change the type of lines plotted and to plot by group rather than sample. We encourage you to explore the documentation to find out more detail.

## Visualizing enrichment in differentially enriched regions

Previously, we had evaluated differential enrichment between the two factors in our study. We had found **almost all of the peaks that were identfied were specific to Nanog and only one region that had significantly higher enrichment in Pou5f1**. We can use the BED files we generated with DiffBind as input to `deepTools` and visualize enrichment in those regions to evaluate the differences in read density.

* Open up `FileZilla` and **copy over the BED files to O2** in`~/chipseq/results/visualization`:

<img src="../img/filezilla_diffbind.png">

Now we can use some of the `deepTools` commands we had explored previously. **Note that we have changed the command from `reference-point` to `scale-regions`.** In the `scale-regions` mode, all regions in the BED file are stretched or shrunken to the length in bases indicated by the user (`--regionBodyLength`).

<img src="../img/computeMatrix_modes.png" width="600">

Let's **start with Nanog file which contains 33 regions** that were identified as increased in enrichment compared to Pou5f1. The plot confirms what we had expected, that is, Pou5f1 don't have much read depth in these regions. 

```bash

 $ computeMatrix scale-regions \
-R ~/chipseq/results/visualization/Nanog_enriched.bed \
-S /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Pou5f1*.bw /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Nanog*.bw \
--skipZeros \
-p 6 \
--regionBodyLength 2000 \
-a 500 -b 500 \
-o ~/chipseq/results/visualization/matrixAll_Nanog_binding_sites.gz


$ plotProfile -m visualization/matrixAll_Nanog_binding_sites.gz \
-out visualization/figures/Allsamples_NanogSites_profile.png \
--perGroup  --plotTitle "" \
--samplesLabel "Pou5f1-Rep1" "Pou5f1-Rep2" "Nanog-Rep1" "Nanog-Rep2" \
-T "Nanog only binding sites"  -z "" \
--startLabel "" \
--endLabel "" \
--colors red red darkblue darkblue
```

<img src="../img/Allsamples_NanogSites_profile.png" width="500">

With **Pou5f1, remember we only had one region**. We are still able to plot this data but you will notice that it is a bit more boxy in nature. This is because values are not being averaged over multiple regions.

```bash

 $ computeMatrix scale-regions \
-R ~/chipseq/results/visualization/Pou5f1_enriched.bed \
-S /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Pou5f1*.bw /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Nanog*.bw \
--skipZeros \
-p 6 \
--regionBodyLength 2000 \
-a 500 -b 500 \
-o ~/chipseq/results/visualization/matrixAll_Pou5f1_binding_sites.gz 


$ plotProfile -m visualization/matrixAll_Pou5f1_binding_sites.gz \
-out visualization/figures/Allsamples_Pou5f1Sites_profile.png \
--perGroup  --plotTitle "" \
--samplesLabel "Pou5f1-Rep1" "Pou5f1-Rep2" "Nanog-Rep1" "Nanog-Rep2" \
-T "Pou5f1 only binding sites"  -z "" \
--startLabel "" --endLabel "" \
--colors red red darkblue darkblue
```

<img src="../img/Allsamples_Pou5f1Sites_profile2.png" width="500">





***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
