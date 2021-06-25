---
title: "Experimental design considerations and understanding the ChIP-seq workflow"
author: "Mary Piper, Radhika Khetani, Meeta Mistry"
date: "March 14th, 2018"
---

Approximate time: 45 minutes

## Learning Objectives

- Guidelines for sequencing depth (depending on signal profile), replicates, and selecting input
- Take content from chipseq workflow slide deck (decision tree)

## Introduction to ChIP-seq
Chromatin immunoprecipitation (ChIP) experiments are performed to identify DNA bound to specific (chromatin) proteins of interest. The first step involves isolating the chromatin and immunoprecipitating (IP) fragements with an antibody against the protein of interest. In ChIP-seq, the immunoprecipitated DNA fragments are then sequenced, followed by identification of enriched regions of DNA or peaks. These peak calls can then be used to make biological inferences by determining the associated genomic features and/or over-represented sequence motifs. 

![chipseq_overview](../img/chipseq_overall.png)
