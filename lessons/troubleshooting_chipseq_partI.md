---
title: "Troubleshooting the ChIP-seq workflow: From sequence reads to peak calls"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu, Will Gammerdinger"
date: "September 15th, 2021"
---

Approximate time: 45 minutes

## Learning Objectives
* Understand the quality checks required at each step of the workflow
* Learn how to troubleshoot issues in the workflow


## From Sequence reads to Peak calls
In this workshop we covered the steps in the first half of the ChIP-seq workflow where we go from raw sequence reads through to peak calls. We have discussed each of the steps in detail, outlining the tools involved and the file formats encountered. In this lesson, we revisit the quality checks associated with each step and summarize the main points to take away.

<p align="center">
<img src="../img/chipseq_deeptoolsworkflow_sept2021.png" width="600">
</p>


## Quality control of sequence reads

The quality checks at this stage in the workflow include:
* Checking the **quality of the base calls** to ensure that there were no issues during sequencing. 
    * Sudden drops in quality could indicate a techincal problem with the sequencer, as shown below.
* Evaluating **duplication levels** for sequences, to check that libraries have not been over amplified. 
    * In a normal, diverse library you should see most sequences occurring only once or twice. In the example below, the unique sequences make up only a small proportion of the library, with low level duplication accounting for most of the sequences. Post-filtering, you will be left with very little usable data.
* Exploring the **sequence length distribution**, as variable lengths can indicate adapter contamination.  
    * Adapter contamination arises from fragments of the sequencing library that are shorter than the read length itself. It is also possible for adapters to be ligated with no insert (a consequence of low DNA or the amount of adapters too high), resulting in only adapter to be sequenced.
    

### Examples of bad quality data

<p align="center">
<img src="../img/badplots_fastqc.png" width="700">
</p>

### Good quality data masquerading as bad quality 

Quality checks include looking for modules of the FastQC report which may report as bad quality (for any other NGS data), but indicate good quality ChIP-seq data:
* Expect to see **over-represented sequences** in the immunoprecipitation sample, because that’s exactly what you’re doing - enriching for particular sequences based on binding affinity.
* High levels of **k-mer content** on either end of the reads suggests is some enrichment of K-mers amongst the fragments. Similar to the point above, we are expecting some similarity in the sequences represented amongst our reads.

<p align="center">
<img src="../img/kmer_fastqc_wt1.png" width="600">
</p>

## Alignment quality

The quality checks at this stage in the workflow include:
* Checking the total percent of reads aligning to the genome 
* Determine the percent of duplicate reads 
* Determining the percent uniquely mapping reads
* Identify percent of reads mapping in blacklist regions
* Checking percent of paired-end reads that are properly paired

Wt-sample 1: low mapping rate:

```
42461605 reads; of these:
  42461605 (100.00%) were unpaired; of these:
    21595346 (50.86%) aligned 0 times
    9875310 (23.26%) aligned exactly 1 time
    10990949 (25.88%) aligned >1 times
49.14% overall alignment rate
```
What went wrong? FASTQC report looks fine

* Try to take the unmapped reads and blast them? Look whether it's a high level of contamination from another organism
* mappability problem? Does bowtie2 have a limit to the number of multimappers i.e > some number considered 'not aligned'
* 

## Peak quality /ChIP quality

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
