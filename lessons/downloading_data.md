---
title: "Downloading ChIP-seq data from external users"
author: "Meeta Mistry, Will Gammerdinger"
date: "September 19th, 2023"
---

Approximate time: 15 minutes

## Learning Objectives
* Download ChIP-seq and reference sequence data used in this workshop for external users

## SRA Toolkit for the Raw Data

To download the files used in this experiment, you will need to use the SRA Toolkit. Once installed, you can use the following commands:

```
cd chipseq_workshop/raw_data/

# For WT Sample 1 ChIP
fastq-dump --accession SRR6823762
fastq-dump --accession SRR6823763
cat SRR6823762.fastq SRR6823763.fastq | gzip > wt_sample1_chip.fastq.gz
rm SRR6823762.fastq SRR6823763.fastq

# For WT Sample 1 Input
fastq-dump --accession SRR6823764
fastq-dump --accession SRR6823765
cat SRR6823764.fastq SRR6823765.fastq | gzip > wt_sample1_input.fastq.gz
rm SRR6823764.fastq SRR6823765.fastq

# For WT Sample 2 ChIP
fastq-dump --accession SRR6823766
fastq-dump --accession SRR6823767
cat SRR6823766.fastq SRR6823767.fastq | gzip > wt_sample2_chip.fastq.gz
rm SRR6823766.fastq SRR6823767.fastq

# For WT Sample 2 Input
fastq-dump --accession SRR6823768
fastq-dump --accession SRR6823769
cat SRR6823768.fastq SRR6823769.fastq | gzip > wt_sample2_input.fastq.gz
rm SRR6823768.fastq SRR6823769.fastq

# For KO Sample 1 ChIP
fastq-dump --accession SRR6823770
cat SRR6823770.fastq | gzip > ko_sample1_chip.fastq.gz
rm SRR6823770.fastq

# For KO Sample 1 Input
fastq-dump --accession SRR6823771
cat SRR6823771.fastq | gzip > ko_sample1_input.fastq.gz
rm SRR6823771.fastq

# For KO Sample 2 ChIP                                                                    
fastq-dump --accession SRR6823772
cat SRR6823772.fastq | gzip > ko_sample2_chip.fastq.gz
rm SRR6823772.fastq

# For KO Sample 2 Input
fastq-dump --accession SRR6823773
cat SRR6823773.fastq | gzip > ko_sample2_input.fastq.gz
rm SRR6823773.fastq
```

## Setting up the reference genome

Let's make a directory to hold our reference genome and move inside it:

```
mkdir chipseq_workshop/reference_genome/
cd chipseq_workshop/reference_genome/
```

You can download the *Mus musculus* reference genome from iGenomes:

```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
```

Unpack the tarball with:

```
tar -xvzf Mus_musculus_UCSC_mm10.tar.gz
```

The iGenomes bundle comes with a pre-indexed genome that is here:

```
ls Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa
```

When we do the genome alignment step with Bowtie2, use this reference instead of the reference that the materials use.

Now, you should be all set to carry out the workflow!
