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


## Experimental design considerations for ChIP-seq

Demonstrate how to design an RNA-seq experiment that avoids confounding and batch effects


Chromatin immunoprecipitation (ChIP) experiments are performed to identify DNA bound to specific (chromatin) proteins of interest. The first step involves isolating the chromatin and immunoprecipitating (IP) fragements with an antibody against the protein of interest. In ChIP-seq, the immunoprecipitated DNA fragments are then sequenced, followed by identification of enriched regions of DNA or peaks. These peak calls can then be used to make biological inferences by determining the associated genomic features and/or over-represented sequence motifs. 

![chipseq_overview](../img/chipseq_overall.png)
