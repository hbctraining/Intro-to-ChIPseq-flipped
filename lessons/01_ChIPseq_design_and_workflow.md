---
title: "Experimental design considerations and understanding the ChIP-seq workflow"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Shannan Ho Sui"
date: "September 15th, 2021"
---

Approximate time: 45 minutes

## Learning Objectives

* Review experimental design considerations for a ChIP experiment
* Understand the steps involved in a ChIP-seq analysis workflow


## Considerations for ChIP-seq

In chromatin immunoprecipitation (ChIP) experiments, a transcription factor, cofactor, or other chromatin protein of interest is enriched by immunoprecipitation from cross-linked cells, along with its associated DNA. In ChIP-seq, the immunoprecipitated DNA fragments are then sequenced, followed by identification of enriched regions of DNA or peaks. These peak calls can then be used to make biological inferences by determining the associated genomic features and/or over-represented sequence motifs. 

<p align="center">
<img src="../img/chipseq_overall.png" width="400">
</p>

ChIP-seq has now been widely used for many transcription factors, histone modifications, chromatin modifying complexes, and other chromatin-associated proteins in a wide variety of organisms. As such, **there is much diversity in the way ChIP-seq experiments are designed and the way analyses are executed.**

In this lesson, we describe a few simple guidelines for setting up a ChIP-seq experiment and outline the analysis workflow.

## Experimental design considerations

When starting out with your experiment, there are many things to think about. We have highlighted some of the important points in the previous lecture and within this lesson, but we also encourage you to peruse the [ENCODE guidelines and practices for ChIP-seq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). Although it was published in 2012, much of the information is still very valid and used in practice today.

> #### What is ENCODE?
> The ENCODE (Encyclopedia of DNA Elements) Project was planned as a follow-up to the Human Genome Project. It aims to identify all functional elements in the human genome. Coinciding with the completion of the Human Genome Project in 2003, the ENCODE Project began as a worldwide effort involving more than 30 research groups and more than 400 scientists.
> 
> Over the years, the ENCODE Consortium has become involved with additional projects whose goals run in parallel. A popular one is the modENCODE (MODel organism ENCyclopedia Of DNA Elements) project, targeting the identification of functional elements in selected model organism genomes, specifically Drosophila melanogaster and Caenorhabditis elegans.

Below, is the **decision tree** that was presented earlier in the workshop (lecture). From this, we highlight some important points below:

<p align="center">
<img src="../img/expt_decisiontree.png" width="900">
</p>


### Starting material
Ensure that you have a sufficient amount of starting material because the ChIP will only enrich for a small proportion. For a standard protocol, you want approximately 2 x 10^6 cells per immunoprecipitation. If it is difficult to obtain that many samples from your experiment, consider using low input methods. Ultimately, higher amounts of starting material yield more consistent and reproducible protein-DNA enrichments.

> NOTE on pooling to increase the amount of starting material. Is it okay? Under what circumstances?


### Quality control of your ChIP

Your ChIP experiment is only as good as your antibody! The more specific the antibody, the more robust and accurate your results will be. Antibody deficiencies are of two main types: poor reactivity against the intended target and/or cross-reactivity with other DNA-associated proteins. Numerous antibodies have been shown to work in ChIP; nevertheless, **it is best to test the antibody with the specific set of cells that you are working with**. Here, we boil it down to the following key points:

* Test your antibody with the use of a **Western blot**. These are performed on protein lysates from either whole-cell extracts, nuclear extracts, chromatin preparations, or immunoprecipitated material. 
* Check a few regions by **qPCR to confirm that the enrichment worked**. This is performed on the immunoprecipitated material, before sending it for sequencing. 
    * You can also check a region of DNA that you do not expect to be enriched and thus do not expect to be amplified by qPCR, to show that your ChIP is specific (negative control)
* If you don't have any known targets for your protein, run a **postive control IP**. Histone H3 or H3k4me3 usually work very well. Since there is loads of H3K4me3 present at most TSSs you could design primers against the promoter of a housekeeping gene. If you have a good signal present, you will at least know the protocol is working well.

> **NOTE:** The authors of this study also included a positive control sample using an antibody against p300 to test the protocol, although the data is not included here. The p300 protein [has been shown](https://pubmed.ncbi.nlm.nih.gov/19212405/) to have binding sites in the cortex. 

### Input control

A ChIP-Seq peak should be compared with the same region of the genome in a matched control sample because only a fraction of the DNA in our ChIP sample corresponds to actual signal amidst background noise. 

There are a number of **artefacts that tend to generate pileups of reads that could be interpreted as a false positive peaks**. These include:
* Open chromatin regions that are fragmented more easily than closed regions due to the accessibility of the DNA
* The presence of repetitive sequences
* An uneven distribution of sequence reads across the genome due to DNA composition
* ‘hyper-ChIPable’ regions: loci that are commonly enriched in ChIP datasets. Certain genomic regions are more susceptible to immunoprecipitation, therefore show increased ChIP signals for unrelated DNA-binding and chromatin-binding proteins.

There are two kinds of controls that can be used for ChIP-seq: **IgG control** and **input control**. 

**IgG control** is DNA resulting from a “mock” ChIP with Immunoglobulin G (IgG) antibody, which binds to non-nuclear antigen. IgG is important for identifying non-specific binding of the beads used for the immunoprecipitation. However, if too little DNA is recovered after immunoprecipitation, the sequencing library will be of low complexity and binding sites identified using this control could be biased. 

**Input control** is DNA purified from cells that are cross-linked, and fragmented, but without adding any antibody for enrichment.
Input samples can account for variations in the fragmentation step of the ChIP protocol as certain regions of the genome are more likely to shear than others based upon their structure and GC content. Input controls are **more widely used** to normalize signals from ChIP enrichment.

> #### Do I need to have one input sample for each IP sample in my dataset?
> The ENCODE guidelines [[Landt, et al, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/)] states _"If cost constraints allow, a control library should be prepared from every chromatin preparation and sonication batch, although some circumstances can justify fewer control libraries. Importantly, a new control is always performed if the culture conditions, treatments, chromatin shearing protocol, or instrumentation is significantly modified."_
> 
> The short answer is, yes having an input sample matched with each IP is ideal. However, studies have shown input replicates to have strong reproducibility in some cases. Thus, if there are constraints with budget or obtaining enough sample, having one input per sample group can suffice. 


### Replicates

As with any high-throughput experiment, a single assay is often subject to a substantial amount of variability. Thus, it is highly recommended to setup your experimental design with a **minimum of 2-3 biological replicates**. Presumably, two replicates measuring the same underlying biology should have high consistency but that is not always the case. Having replicates allow you to evaluate concordance of peaks and identify a set of reproducible enriched regions with greater confidence. If you have multiple sample groups and are planning a differential enrichment analysis, increasing the number of replicates will give you more statistical power to find changes between groups.

> #### Do we see batch effects in ChIP-seq data?
> Typically, batch effects are not as big of a concern with ChIP-seq data. However, it is best to run everything in parallel as much as possible. If you only have a single sample group, it should be more feasible to prepare all samples together (since there are fewer). For multiple sample groups, if you are not able to process all samples together, split replicates of the different sample groups across batches. This way you avoid any potential confounding.


### Sequencing considerations

Below we list some general guidelines and things to think about when sending your samples to the sequencing facility:

#### Read length
* Read length shoudl be between 50- to 150-bp
* Single-end reads is sufficient in most cases
    * Longer reads and paired-end reads will improve mappability
    * Paired-end is good (and necessary) for allele-specific chromatin events, and investigations of transposable elements
* Balance cost with value of more informative reads
    * i.e. spending money on replicates is more important than longer reads or paired-end

#### Sequencing depth

* Narrow peak profiles
     * **Mammalian cells**; ENCODE suggests a minimum of 10 million uniquely mapped reads. For standard transcription factors we recommend between **20-40 million total read depth**
     * **Worms and flies**; ENCODE suggests a minimum of **2 million uniquely mapped reads**. We recommend between **4-8 million total read depth**.
     
* Broad peak profiles
     * Generally require a higher sequence depth
     * **Mammalian cells** require a **minimum of 40M total read depth; more is better** for detecting some histone marks
     * **Worms and flies; less is known** and so the numbers vary across studies. We suggest a **minimum of 8M total read depth.**

* Sequence the input controls to equal or higher depth than your ChIP samples

### Resources
* [Thermofisher Step-by-step guide to a successfuk ChIP experiment](https://www.thermofisher.com/us/en/home/life-science/antibodies/antibodies-learning-center/antibodies-resource-library/antibody-application-notes/step-by-step-guide-successful-chip-assays.html)
* [Profiling of transcription factor binding events by ChIP-seq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5544034/)

## Understanding the ChIP-seq analysis workflow

In the same way that the experimental design setup is a process, the analysis also takes some thought and decision making. Below we have depicted and end-to-end analysis workflow for ChIP-sequencing data. **Our focus for this workshop will be on the top half of this workflow diagram**.

All of the software required to get us from **raw sequence reads to peak calls** are command-line tools and accessible on O2, the HMS-RC high performance cluster. As we encounter each of the tools, we will describe it's purpose and how it functions. Where applicable we can delve a bit deeper to understand the inner workings of the algorithm, such that we can thororughly understand the output. Each step of the workflow will require a specific file format. You will notice that we have annotated these file formats on the workflow below. We will describe these file formats in more detail in the respective lessons.

<p align="center">
<img src="../img/chipseq_fullworkflow_sept2021.png" width="700">
</p>

> **NOTE:** Boxes in green cover topics that will not be covered here, but will be a focus of the ChIP-seq Part II workshop (_material development in progress_). It is also worth noting that the green boxes also represent analysis steps which require a working knowledge of R. The color fade indicates that there are some parts of "Peak Call QC" which are covered in this workshop with the use of command-line tools, and other methods which are R-based and will be covered in Part II. 

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*



