---
title: "Experimental design considerations and understanding the ChIP-seq workflow"
author: "Mary Piper, Radhika Khetani, Meeta Mistry"
date: "March 14th, 2018"
---

Approximate time: 45 minutes

## Learning Objectives

* Provide guidelines for sequencing depth
* Explain the protocol for creating an appropriate input
* Describe the importance of replicates for ChIP-seq experiments
* Understand the possible routes of analysis in the ChIP-seq workflow


## Considerations for ChIP-seq

In chromatin immunoprecipitation (ChIP) experiments, a transcription factor, cofactor, or other chromatin protein of interest is enriched by immunoprecipitation from cross-linked cells, along with its associated DNA. In ChIP-seq, the immunoprecipitated DNA fragments are then sequenced, followed by identification of enriched regions of DNA or peaks. These peak calls can then be used to make biological inferences by determining the associated genomic features and/or over-represented sequence motifs. 

<p align="center">
<img src="../img/chipseq_overall.png" width="400">
</p>

ChIP-seq has now been widely used for many transcription factors, histone modifications, chromatin modifying complexes, and other chromatin-associated proteins in a wide variety of organisms. As such, there is much diversity in the way ChIP-seq experiments are designed and the way analyses are executed. In this lesson, we describe a few simple guidelines for setting up a ChIP-seq experiment and where applicable discriminate between requirements for narrow versus broad peaks.

