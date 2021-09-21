---
title: "Peak calling with MACS2"
author: "Meeta Mistry, Jihe Liu, Radhika Khetani, Mary Piper, Will Gammerdinger"
date: "September 20th, 2021"
---

Contributors: Meeta Mistry, Jihe Liu, Radhika Khetani, Mary Piper, Will Gammerdinger

Approximate time: 60 minutes

## Learning Objectives

* Describe the different components of the MACS2 peak calling algorithm
* Describe the parameters involved in running MACS2
* List and describe the output files from MACS2

## Peak Calling

Peak calling, the next step in our workflow, is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment. 

<p align="center">
<img src="../img/chipseq_peakworkflow_sept2021.png" width="700">
</p>

For ChIP-seq experiments, what we observe from the alignment files is a **strand asymmetry with read densities on the +/- strand, centered around the binding site**. The 5' ends of the selected fragments will form groups on the positive- and negative-strand. The distributions of these groups are then assessed using statistical measures and compared against background (input or IgG samples) to determine if the site of enrichment is likely to be a real binding site.

<p align="center">
<img src="../img/plos_chipseq_arrow.png" width = "700">
</p>

*Image source: [Wilbanks and Faccioti, PLoS One 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011471)*

There are various tools that are available for peak calling. Peak calling algorithms are often specialized in identifying one of **two types of enrichment**: broad peaks or narrow peaks. There are also many tools out there that are capable of handling both types of profiles, and have specific methods for each. As such, it is good to **have some idea about what type of binding profile you are expecting when choosing your peak caller and/or the specific methods** to run. For more detail on the different types of binding profiles, please refer to the discussion from [an earlier lesson](01_ChIPseq_design_and_workflow.md#types-of-binding-profiles).

> #### What if I am not sure what type of binding profile to expect?
> In this scenario, we advise peak calling for both narrow and broad profiles. This will require additional time for some exploration to determine what is best for your data.


## MACS2
A commonly used tool for identifying binding sites is named [Model-based Analysis of ChIP-seq (MACS)](https://github.com/taoliu/MACS). The [MACS algorithm](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. Although it was developed for the detection of transcription factor binding sites (narrow peaks), it is also suited for larger regions (broad peaks). 

We will be using the newest version of this tool, **MACS2**. The **underlying algorithm for peak calling remains the same as before**, but it comes with some enhancements in functionality. The MACS/MACS2 workflow is depicted below. In this lesson, we will describe the steps in more detail.

> **NOTE:** While MACS can be used to call peaks without an input control, we advise against this. The control sample increases specificity of the peak calls, and without it you will find many false positive peaks identified. 

<p align="center">
<img src="../img/macs_workflow.png" width="400">
</p>
	
### Removing redundancy

MACS provides different options for dealing with **duplicate tags** at the exact same location, that is tags with **the same coordination and the same strand**. The default is to keep a single read at each location. The `auto` option, which is very commonly used, tells MACS to calculate the maximum tags at the exact same location based on binomal distribution using 1e-5 as the pvalue cutoff. An alternative is to set the `all` option, which keeps every tag. If an `integer` is specified, then at most that many tags will be kept at the same location. This redundancy is consistently applied for both the ChIP and input samples.

> _We do not need to worry about this option, since **we filtered out the duplicates during the [post-alignment filtering step](05_filtering_BAM_files.md).**_


### Modeling the shift size

The tag density around a true binding site should show a **bimodal enrichment pattern** (or paired peaks). MACS takes advantage of this bimodal pattern to empirically model the shifting size, thus better locating the precise binding sites.

To identify the shift size:

1. MACS scans the whole sample **searching for all highly significant enriched regions**. *This is done only using the ChIP sample!* 
   * These regions are identified by MACS sliding two `bandwidth` windows across the genome to find regions with **tags more than `mfold` enriched relative to a random tag genome distribution**. 
2. MACS randomly **samples 1,000 of these high-quality peaks** identified in #1. 
3. For these 1,000 peaks, MACS separates their positive and negative strand tags and aligns them by the midpoint between their centers. The **distance between the modes of the two peaks in the alignment is defined as 'd'** and represents the estimated fragment length. 
4. MACS **shifts all reads in the sample by d/2** toward the 3' ends to the most likely protein-DNA interaction sites

<p align="center">
<img src="../img/peak_shift3.png" width="400">
</p>

### Scaling libraries

For experiments in which sequence depth differs between input and treatment samples, MACS linearly scales the **total control tag count to be the same as the total ChIP tag count**. The default behaviour is for the larger sample to be scaled down.

### Effective genome length

To calculate λBG (a parameter discussed in "Peak detection" below) from tag count, MACS requires the **effective genome size** or the size of the genome that is mappable. Mappability is related to the uniqueness of the k-mers at a  particular position the genome. Low-complexity and repetitive regions have low uniqueness, which means low mappability. Therefore we need to provide the effective genome length to **correct for the loss of true signals in low-mappable regions**.

<p align="center">
<img src="../img/mappable.png" width="300">
</p>

> #### How do I obtain the effective genome length?
> The MACS software has some pre-computed values for commonly used organisms (human, mouse, worm and fly). If you wanted you could compute a more accurate values based on your organism and build. The [deepTools docs](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) has additional pre-computed values for more recent builds but also has some good materials on how to go about computing it.

### Peak detection

After MACS shifts every tag by *d/2*, it then slides across the genome using a window size of *2d* to find candidate peaks. The tag distribution along the genome can be modeled by a Poisson distribution. The Poisson is a one parameter model, where the parameter **λ is the expected number of reads in that window**.

<p align="center">
<img src="../img/peak_detection.png" width="300">
</p>

Instead of using a uniform λ estimated from the whole genome (computed using only the input control), MACS uses a λlocal defined for each candidate peak. The λlocal parameter is deduced by **taking the maximum λ value across the λ computed across various window sizes** (as shown below): 

<p align="center">
<img src="../img/lambda.png" width="300">
</p>

**λlocal = max(λ300bp, λ1KB, λ5KB, λ10KB, λBG).** 

> λBG represents the background λ estimated using the whole genome (i.e. the largest window size)

In this way, lambda captures the influence of local biases, and is **robust against occasional low tag counts at small local regions**. Possible sources for these biases include local chromatin structure, DNA amplification and sequencing bias, and genome copy number variation.

Next, a Poisson distribution p-value is computed based on λ. A region is considered to have a significant tag enrichment if the p-value < 1e-5. Any overlapping enriched peaks are merged into a single peak.

### Estimation of false discovery rate

Each peak is considered an independent test. Therefore, when we encounter thousands of significant peaks detected in a sample, we have a multiple testing problem. In MACSv1.4, the FDR was determined empirically by exchanging the ChIP and control samples. However, in MACS2, p-values are now corrected for multiple comparison using the **Benjamini-Hochberg correction**.

## Running MACS2 

To run MACS2, we will first load the macs2 module along with any dependencies:

```bash
$ module load gcc/6.2.0  python/2.7.12 macs2/2.1.1.20160309
```

We will also need to create a directory for the output generated from MACS2:

```bash
# Create macs2 directory in results
$ mkdir -p ~/chipseq_workshop/results/macs2
```

Now change directories to the `results` folder:

```bash
$ cd ~/chipseq_workshop/results/
```

Since we only created a filtered BAM file for a single sample, **we will use the BAM files we have generated for you**. Rather than copying them over, we will have you point to them within your peak calling command.


### MACS2 parameters

There are seven [major functions](https://github.com/taoliu/MACS#usage-of-macs2) available in MACS2 serving as sub-commands. We will only cover `callpeak` in this lesson, but you can use `macs2 COMMAND -h` to find out more, if you are interested.

`callpeak` is the main function in MACS2 and can be invoked by typing `macs2 callpeak`. If you type this command without parameters, you will see a full description of command-line options. Here is a short list of the commonly used ones: 

**Input file options**

* `-t`: The ChIP data file (this is the only REQUIRED parameter for MACS)
* `-c`: The control or mock data file
* `-f`: format of input file; Default is "AUTO", which will allow MACS to decide the format automatically.
* `-g`: mappable genome size, which is defined as the genome size that can be sequenced; some precompiled values provided.

**Output arguments**

* `--outdir`: MACS2 will save all output files into speficied folder for this option
* `-n`: The prefix string for output files
* `-B/--bdg`: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

**Shifting model arguments**

* `-s`: size of sequencing tags. By default, MACS will use the first 10 sequences from your input treatment file to determine it
* `--bw`: The bandwidth, which is used to scan the genome ONLY for model building. Can be set to the expected sonication fragment size.
* `--mfold`: upper and lower limit for model building

**Peak calling arguments**

* `-q`: q-value (minimum FDR) cutoff
* `-p`: p-value cutoff (instead of q-value cutoff)
* `--nolambda`: do not consider the local bias/lambda at peak candidate regions
* `--broad`: broad peak calling

> **NOTE:** Relaxing the q-value does not behave as expected in this case, since it is partially tied to peak widths. Ideally, if you relaxed the thresholds, you would simply get more peaks. But with MACS2, relaxing thresholds also results in wider peaks.

Now that we have a feel for the different ways we can tweak our command, let's set up the command for each of our wildtype replicates:

```
macs2 callpeak -t /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_chip_final.bam \
    -c /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_input_final.bam \
    -f BAM -g mm \
    -n wt_sample1 \
    --outdir macs2 2> macs2/wt_sample1_macs2.log

$ macs2 callpeak -t /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip_final.bam \
    -c /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_input_final.bam \
    -f BAM -g mm \
    -n wt_sample2 \
    --outdir macs2 2> macs2/wt_sample2_macs2.log
```

The tool is quite verbose, and normally you would see lines of text being printed to the terminal, describing each step that is being carried out. We have captured that information into a log file using `2>` which re-directs the stadard error into a file. **You can use `less` to look at the log file and see what information is being reported.**
	
Move the log files to the `log` directory we had created during our project setup:

```bash
$ mv macs2/*.log ../logs/
```

## MACS2 Output files

Change directories into `macs2`, and list the output files that we have generated.

```bash
$ cd macs2/
	
$ ls -lh
```

There should be 4 output files in the results directory for each sample (2 replicates), so a total of 8 files:

* `_peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
* `_peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment (the ratio between the ChIP-seq tag count and λlocal)
* `_summits.bed`: The location in the peak with the highest fragment pileup. These are the predicted precise binding location and recommended to use for motif finding.
* `_model.R`: an R script which you can use to produce a PDF image about the model based on your data and cross-correlation plot


Let's first obtain a summary of how many peaks were called in each sample. We can do this by counting the lines in the `.narrowPeak` files:

```bash
$ wc -l *.narrowPeak
```
	
*** 

**Exercise:**

1. Using the BAM files listed below, use MACS2 to call peaks on the KO samples.

```bash
# KO ChIP BAM files
/n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/ko_sample1_chip_final.bam
/n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/ko_sample2_chip_final.bam

# KO input BAM files
/n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/ko_sample1_input_final.bam
/n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/ko_sample2_input_final.bam
```

2. How many peaks do you get for each sample? How does this compare with the WT?

***

In the next lesson, we will delve deeper into the output files and gain an understanding of the different file formats.

### Other peak calling software
We have used this lesson to describe to you the inner workings of the MACS2 peak caller. This is a very popular tool in the field, however there are many other peak callers and depending on the data you are working with, it can be worth exploring more. Below, we provide the names of some other peak calling software, but this is by no means an exhaustive list.

* [SPP](https://www.encodeproject.org/software/spp/): an R package, that is implemented in the ENCODE processing pipeline. Best for narrow peak calling. 
* [EPIC](https://github.com/biocore-ntnu/epic2): ideal for broad peak calling (a re-implentation of an older tool called SICER)
* [haystack bio](https://github.com/pinellolab/haystack_bio): Epigenetic Variability and Motif Analysis Pipeline


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
