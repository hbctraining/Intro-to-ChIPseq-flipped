---
title: "ATAC-seq data analysis"
author: "Meeta Mistry"
date: "July 7th,  2022"
---

## Assaying open regions in genome 

If we unwind the DNA from the chromosome, we will uncover that a complex of macromolecules exist in an intricately organized manner. The most **basic repeating unit of chromatin is the nucleosome**, consisting of ∼147 bp of DNA wrapped around an octamer of histone proteins.

* There are sections of highly condensed **heterochromatin**, where nucleosomes are packed into chromatin fiber. This is typically regions that are gene-poor and/or transcriptionally silent.

* On the other hand, there is **euchromatin** which is is more loosely wrapped chromatin. These are **regions of chromatin that are more open**, allowing recruitment of RNA polymerase complexes and gene regulatory proteins for active transcription


<p align="center">
<img src="img/chromatin.jpeg" width="500">
</p>

_Image source: ["Creative Diagnostics Blog"](https://www.creative-diagnostics.com/blog/index.php/the-structure-and-function-of-chromatin/)_

## ATAC-seq

A popular approach used to identify open regions of the genome is the **A**ssay for **T**ransposase-**A**ccessible **C**hromatin (ATAC) followed by high throughput sequencing.  The ATAC-Seq method was first [published in 2013](https://www.ncbi.nlm.nih.gov/pubmed/24097267) in the journal Nature Methods by lead researcher Jason Buenrostro in the labs of Howard Chang and William Greenleaf at Stanford University.

### How does it work?

* Utilizes hyperactive Tn5 transposase to insert sequencing adapters into the open chromatin regions 
* Tn5 tagmentation simultaneously fragments the genome and tags the resulting DNA with sequencing adapters
* Amplify and sequence

<p align="center">
<img src="img/atacseq_schematic.png" width="500">
</p>

_Image source: [Buenrostro et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/24097267)_

The method relies on the hyperactive Tn5 transposase that was already being used for tagmentation-based NGS library preparation methods. The authors hypothesized that if a similar approach was used in vivo, the **addition of the adapters would mainly take place in open chromatin regions, where no steric hindrance of the transposase would occur**, allowing the enzyme to preferentially access these regions.

### The result

ATAC-seq simultaneously assesses three chromatin mechanisms together in one assay:

* Maps open chromatin
* Transcription factor occupancy
* Nucleosome occupancy


