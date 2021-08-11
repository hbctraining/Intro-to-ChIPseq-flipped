---
title: "Visualization of peaks"
author: "Meeta Mistry, Jihe Liu"
date: "Aug 11th, 2021"
---

Approximate time: 80 minutes

**Link to issue describing the modifications to be made:** https://github.com/hbctraining/Intro-to-ChIPseq-flipped/issues/12

## Learning Objectives

* Visualizing enrichment patterns at particular locations in the genome

## Profile plots and heatmaps

After creating the bigWig files, we are ready to perform data visualization. We first need to prepare an intermediate file that will be used with the `plotHeatmap` and `plotProfile` commands.

<p align="center">
<img src="../img/computeMatrix_overview.png" width="700">
</p>

The `computeMatrix` command accepts multiple bigWig files and multiple region files (BED format) to create a count matrix -- the intermediate file. The command can also filter and sort regions according to their scores. The region file will be the BED file we generated for the final peaks, and the bigWig files will be those generated from the last lesson. Additionally, we specify a window of +/- 4000 bp around the reference point of genes (`-b` and `-a`). The [reference point](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#Optional%20arguments_repeat1) for the plotting could be either the region start (TSS), the region end (TES) or the center of the region. Here, we use the center of the region (the default is TSS). For each window, `computeMatrix` will calculate scores based on the read density values in the bigWig files.

Let's create a matrix for the wt sample:

```bash

computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S ~/chipseq_workshop/results/visualization/bigWig/wt_sample1_chip.bw ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-p 6

```

> **NOTE:** Typically, the genome regions are genes, and can be obtained from the [UCSC table browser](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables). Alternatively, you could look at other regions of interest that are not genomic feature related (i.e. binding regions from another protein of interest).

With the computed matrix, we could now create a **profile plot**, which is essentially a density plot that evaluates read density across all reference points. For the wt samples, we can see that **sample2 has higher amount of signal at the reference point compared to sample1**. 

```bash

# Create figures directory under visualization
mkdir ~/chipseq_workshop/results/visualization/figures

# Plot the profiles
plotProfile -m ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot1_wt.png \
--regionsLabel "" \
--perGroup \
--colors red blue \
--samplesLabel "PRDM16_sample1" "PRDM16_sample2" \
--refPointLabel "PRDM16 binding sites"

```

<p align="center">
<img src="../img/09_plot1_wt.png" width="500">
</p>


> **NOTE:** Both `plotProfile` and `plotHeatmap` have many options, including the ability to change the type of lines plotted, and to plot by group rather than sample. We encourage you to explore the documentation to find out more detail.

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
