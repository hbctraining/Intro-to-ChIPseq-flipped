## Experimental design consideration

1. Which of the following is False:
  - **To normalize signals, IgG control is used more often than input control in a ChIPseq analysis**
  - A minimum of 2-3 biological replicates are required for a ChIPseq experiment
  - Techniques to validate antibody efficiency and specificity include western blot and qPCR

1. Which of the following is True about the peaks:
  - CTCF and RNA Polymerase are examples that present broad peaks
  - Narrow peaks are often associated with histone modifications
  - **For a same organism, broad peaks require more sequencing depths, compared to narrow peaks**

## Bowtie2 alignment

1. True or False: for our dataset, we use the local alignment mode, instead of global, when performing Bowtie2 alignment.
  - **True**
  - False

1. Which of the following about the BAM is False:
  - It contains the same information as the SAM format
  - **It is the default output of Bowtie2 aligner**
  - It is a binary file, which is not human readable
  - You could convert a file in SAM format to BAM format using "samtools view" command

## Handling peak files

1. Which of the following statement is False:
  - Compared to BED file, narrowPeak file provides additional statistical information
  - BED file always has these three columns: chromosome, start position, end position
  - **Blacklist regions should only be removed after the peak calling**
  - It is best practise to remove peaks that overlap with blacklist regions

1. Which of the following command is used to find overlap peaks between replicates:
  - bedtools merge
  - bedtools genomecov
  - **bedtools intersect**
  - bedtools subtract
