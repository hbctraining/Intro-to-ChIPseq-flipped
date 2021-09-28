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

* Checking the **total percent of reads aligning to the genome**. A mapping rate of 70% or higher is great, and lower than 50% is concerning especially if total read depth is low.
   * Low mapping rates can be a result of low quality reads, contaminating sequences, inappropriate alignment parameters chosen, or a poor quality reference genome. 

> **If my mapping rate is low, do I discard my sample?** Do not discard your sample, rather you will want to:
> 1. Flag the sample as low quality. Keep an eye out for QC metrics later in the workflow for that same sample.
> 2. Troubleshoot the sample. Take the unmapped reads and BLAST the sequences. It's possible you might identify a high level of contamination from another organism.
 
* Determining the **percent uniquely mapping reads**, ideally this would be > 60% of the total read depth of your sample. The higher the percentage the more usable data you have to work with.
   *  High duplication can be a result of over-amplification. 
   *  High number of multi-mapping reads can be due to mappability.

<p align="center">
<img src="../img/mapping_reads.png" width="600">
</p>

_Image source: [Land et, al, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/)_

* Identify **percent of reads mapping in blacklist regions**. This shouldn't be higher than 10% of your high quality mapped (and filtered) read count.

> **NOTE:** For paired-end reads you will also want to checking percent that are properly paired. By default, Bowtie 2 searches for both concordant and discordant alignments, though searching for discordant alignments can be disabled with the `--no-discordant option`.


## Peak quality




***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
