## Experimental design consideration

1. Which of the following statements is FALSE when designing your ChIP-seq experiment:
  - If you have 5ng of IP DNA or less, your may want to consider using a low input ChIP method 
  - **To normalize signals in ChIP-seq analysis, IgG control is used more often than input control**
  - A minimum of 2-3 biological replicates are required for a ChIPseq experiment
  - Techniques to validate antibody efficiency and specificity include western blot and qPCR

1. Which of the following statements about binding site profiles is TRUE:
  - CTCF and RNA Polymerase are examples that present broad peaks
  - Narrow peaks are often associated with histone modifications
  - **For the same organism, broad peaks require more sequencing depths, compared to narrow peaks**
  - You will only find one of two types of binding profiles: narrow or broad
  - All histone marks wil generate broad peaks

## First steps of the workflow: QC and alignment

1. Which of the following FASTQC plots may report as "FAILURE", but is actually a sign of good signal in your ChIP-seq data?
  - Sequence duplication
  - Sequence length distribution 
  - **K-mer content**
  - Per base sequence quality
  - Per base sequence content

1. True or False: For our dataset, we use the "local" alignment mode, instead of "global", when using the Bowtie2 aligner.
  - **True**
  - False

1. Which of the following statements about the BAM file format is FALSE:
  - It contains the same information as the SAM format
  - It is a binary file, which is not human readable
  - You can convert a file in SAM format to BAM format using "samtools view" command
  - **It is the default output of Bowtie2 aligner**

## Handling peak files

1. Which of the following statements is FALSE:
  - **Blacklist regions should only be removed after the peak calling**
  - Compared to a BED file, narrowPeak file provides additional statistical information
  - BED file has at minimum three columns: chromosome, start position, end position
  - It is best practise to remove peaks that overlap with blacklist regions

1. Which of the following command is used to find overlapping peaks between replicates:
  - bedtools merge
  - bedtools genomecov
  - **bedtools intersect**
  - bedtools subtract


## Qualitative assessment of peak enrichment

1. Which of the following file formats is NOT typically used directly for visualizing ChIP-seq data:
  - BED 
  - bedGraph
  - **BAM**
  - bigWig
  - narrowPeak

1. True or False: When evaluating signal profile plots and heatmaps, you do not observe any enrichment around the TSS - this is an indication that your experiment has failed.

1. Which of the following statements is FALSE:
  - A histone mark can provide us information about active promoter regions in the genome
  - The profile plot does not always have to be centered on the TSS of genes
  - **The qualitative assesment allows for defined conclusions about our data.**
  - A heatmap can be used as an alternative visualization for a profile plot, as the same data is being plotted
  - Using public data is helpful when your budget does not allow for additional experiments
