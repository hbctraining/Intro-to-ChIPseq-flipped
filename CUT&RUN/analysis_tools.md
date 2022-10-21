### [CUT&RUNTools](https://github.com/fl-yu/CUT-RUNTools-2.0)
An end-to-end computational pipeline specifically tailored to this technology. Latest version now available for single-cell analysis. CUT&RUNTools is implemented using Python, R, and BASH scripts.

**Workflow/features:**
* Takes paired-end sequencing readFASTQ files as the input and performs a set of analytical steps
* Trimming of adapter sequences (Trimmomatic). A two-step read trimming process to improve the quality (K-seq).
* Alignment to the reference genome (Bowtie2). Turning on dovetail alignment, designed to accept alignments for paired-end reads when there is a large degree of overlap between two mates of each pair.
* Size selection: After alignment, fragments are divided into ≤ 120-bp and > 120-bp fractions. The ≤ 120-bp fraction which is likely to contain TF binding sites.
* Peak calling (MACS2, and now SEACR in 2.0)
* Estimation of cut matrix at single-nucleotide resolution (mostly used for footprinting)
* De novo motif searching (MEME) and motif footprinting analysis, using sequences within 100 bp from the summit of each peak
* Direct binding site identification
* Data visualization


<p align="center">
<img src="img/CRtools_figure1.png" width=500>
</p>


_Image source: ["CUT&RUNTools: a flexible pipeline for CUT&RUN processing and footprint analysis"](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-019-1802-4.pdf)_


### [SEACR](https://seacr.fredhutch.org/)
Peak calling by Sparse Enrichment Analysis for CUT&RUN sequencing data. An analysis strategy that uses the global distribution of background signal to calibrate a simple threshold for peak calling. 

* CUT&RUN data features exceedingly low background and low sequence depth, in comparison with ChIP-seq. 
* ChIP-seq experiments are typically sequenced deeply and thus feature high background, thus most peak calling algorithms designed for this type of data. The sparseness of the background can increase false positives, resulting in any spurious background read being called as a peak.
* Thus, rather than requiring highly sensitive methods to distinguish signal from background noise, **peak calling from CUT&RUN data requires high specificity for true positive peaks**.

<p align="center">
<img src="img/seacr_fig1.png" width=500>
</p>

_Image source: ["Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling"](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)_


**Features of SEACR:**

- model free and empirically data driven and therefore does not require arbitrary selection of parameters from a statistical model
- fast, accurate, scalable and simple to use
- peak calling is based on fragment block aggregation

### [CnRAP (Cut & Run Analysis Pipeline)](https://github.com/mbassalbioinformatics/CnRAP)
An analytical pipeline developed to analyze CUT&RUN data. Inspired by both Henikoff (SEACR) and Orkin (Cut&RunTools) lab pipelines.

