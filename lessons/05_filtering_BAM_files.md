---
title: "Alignment and filtering"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu, Will Gammerdinger"
date: "Aug 10th, 2021"
---

Contributors: Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu, Will Gammerdinger

Approximate time: 45 minutes

## Modifications

**CUT&RUN**

* Duplicate removal is optional in many CUT&RUN analysis approaches (Henikoff, CUT&RUN tools, and STAR protocols). 
  * **Default is to keep duplicates. But, why?** Quoted from CUT&RUNtools paper: "...it is argued that nuclease cleavage of chromatin by its stereotypical nature is influenced by conformation of chromatin and/or nuclease bias, and shorter DNA fragments also increased the likelihood of identical reads that originated from different cells." Basically they say that the nature of the technique will increase likelihood of biological duplicates. Remove duplicates from CUT&RUN experiments with caution. First assess your library complexity and if you do not see an unreasonably high amount of duplication and you know you didn't over-amplify, you might not want to remove.
* No mention of multi-mapper removal in any CUT&RUN analysis approaches either. Maybe because Bowtie2 defaults to search for multiple alignments, yet only reports the best one (i.e. not using the `-k` option)
* Filtering of BAM files by fragment size.
   * After alignment, fragments can be divided into ≤ 120-bp and > 120-bp fractions. For transcription factors or proteins with an expected punctate binding profile, you can use the ≤ 120-bp fraction which is likely to contain binding sites. The range can be increased depending on the protein of interest, and alternatively BAM files without filtering can also be used. 
   * Done using sambamba: 

```bash
sambamba view --format \
  bam --nthreads 6 \
  -F "((template_length > 0 and template_length < 120) or (template_length < 0 and template_length > -120))" $file | samtools view -b > bams_sizeSelect/${s}-sizeSelect.bam
```

**ATAC-seq**

* Uniquely mapping reads is critical. A unique mapping rate over 80% is typical for a successful ATAC-seq experiment. Duplicates and multi-mappers are removed. 
* Filtering mitochondrial reads. The mitochondrial genome, which is more accessible due to the lack of chromatin packaging will result in extremely high read coverage. These reads should be discarded. Since there are no ATAC-seq peaks of interest in the mitochondrial genome, these reads are discarded. The Omni-ATAC method uses detergents to remove mitochondria from the samples prior to sequencing and is another option to deal with this issue.
* Filtering BAM files based on fragment size.
   * Typically, a successful ATAC-seq experiment should generate a fragment size distribution plot with decreasing and periodical peaks corresponding to the nucleosome-free regions (NFR) (< 100 bp) and mono-, di-, and tri-nucleosomes (~ 200, 400, 600 bp, respectively)
   * Fragments from the NFR are expected to be enriched around the transcription start site (TSS). Fragments from nucleosome-bound regions are expected to be depleted at TSS with a slight enrichment of flanking regions around TSS.  Use sambamba code from above to filter out fragments by size. A BAM for NFR, mono-nuc, di-nuc, tr-nuc. Typically the NFR BAM is used for peak calling.  
   * Shift the reads in the BAM file. Reads should be shifted + 4 bp and − 5 bp for positive and negative strand respectively, to account for the 9-bp duplication created by DNA repair of the nick by Tn5 transposase

## Learning Objectives

* Describe the purpose of filtering alignment reads
* Employ filtering on alignment files using sambamba and samtools

## Filtering reads

<p align="center">
 <img src="../img/chipseq_filterworkflow_sept2021.png" width="600">
</p>


A key issue when working with a ChIP-seq data is to **move forward with only the uniquely mapping reads**.  Allowing for multi-mapped reads increases the number of usable reads and the sensitivity of peak detection; however, the number of false positives may also increase [[1]](https://www.ncbi.nlm.nih.gov/pubmed/21779159/). To increase our confidence in peak calling and improve data reproducibility, we need to **filter out both multi-mapping reads and duplicate reads**.

* Multi-mapping reads are reads that are mapping to multiple loci on the reference genome.

<p align="center">
 <img src="../img/Multimapping_reads.png" width="500">
</p>

* Duplicate reads are reads that map at the exact same location, with the same coordinates and the same strand. These duplicates can arise from experimental artifacts, but can also contribute to genuine ChIP-signal.
    * **The bad kind of duplicates:** If initial starting material is low, this can lead to overamplification of this material before sequencing. Any biases in PCR will compound this problem and can lead to artificially enriched regions. 
    * **The good kind of duplicates:** You can expect some biological duplicates with ChIP-seq since you are only sequencing a small part of the genome. This number can increase if your depth of coverage is excessive or if your protein only binds to few sites. If there are a good proportion of biological dupicates, removal can lead to an underestimation of the ChIP signal. 

<p align="center">
 <img src="../img/Duplicate_reads.png" width="500">
</p>

> #### Some additional notes on duplicates
> Most peak calling algorithms also implement methods to deal with duplicate reads. While they are commonly removed prior to peak calling, another option is to leave them now and deal with them later. **Skip the duplicate filtering at this step if**:
> * You are planning on performing a differential binding analysis.
> * You are expecting binding in repetitive regions (also, use paired-end sequencing) 
> * You have included [UMIs](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-molecular-identifiers.html) into your experimental setup.

The older version of Bowtie2 had an argument that allowed us to easily perform filtering during the alignment process. but the latest Bowtie2 does not have this option. As a result, the filtering will be done with the use of a tool called [sambamba](https://lomereiter.github.io/sambamba/). Sambamba is an open source tool that provides methods for working with SAM/BAM files, similar to samtools, except with faster processing times and in some cases added functionality. 

This **lesson will consist of two steps**:

1. Sort BAM files by genomic coordinates (using `samtools`).
2. Filter the reads to keep only uniquely mapping reads (using `sambamba`). This will also remove any unmapped reads.

Before we begin, you will want to make sure you are **logged into O2.** To start an interactive session with 2 cores and 10G of memory (sorting can be memory-intensive) us the command below:

> **NOTE:** Check to see that you are not already on a compute node! Interactive session should always be run from the login node. If the text before your command prompt contains the word 'compute', _do not run the code below_.

```bash
$ srun --pty -p interactive -t 0-2:30 --mem 10G -c 2 --reservation=HBC2 /bin/bash
```

We will also load the required modules for this lesson:

```bash
module load gcc/6.2.0 samtools/1.13 sambamba/0.7.1
```


### 1. Sort BAM files by genomic coordinates

Before we can do the filtering, we need to sort our BAM alignment files by genomic coordinates (instead of by name). To perform the sorting, we could use [Samtools](http://www.htslib.org/), a tool we previously used when coverting our SAM file to a BAM file. 

The command we use this time is `samtools sort` with the parameter `-o`, indicating the path to the output file. The example code looks like below (please do not run):

``` bash
# DO NOT RUN
$ cd ~/chipseq_workshop/results/bowtie2/
$ samtools sort wt_sample2_chip.bam -o wt_sample2_chip_sorted.bam
```

> **NOTE**: To run the above code, you will need the BAM file generated from the [alignment lesson](04_alignment_using_bowtie2.md). If you do not have this file, you could copy over the BAM file to your directory:
>
> ```bash
> $ cp /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip.bam ~/chipseq_workshop/results/bowtie2/wt_sample2_chip.bam
> ``` 

The `samtools sort` code above takes about 7 min to finish. Instead of running it in the class, we have generated the output BAM file. Please copy over the BAM file to your directory:

```bash
$ cd ~/chipseq_workshop/results/bowtie2/
$ cp /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip_sorted.bam ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_sorted.bam
```

We could take a glimpse of the sorted BAM file using `samtools view`:

```bash
$ samtools view wt_sample2_chip_sorted.bam | less
```

### 2. Filter the reads to keep only uniquely mapping reads

Next, we can filter the sorted BAM files to keep only uniquely mapping reads. We will use the `sambamba view` command with the following parameters:

* `-t`: number of threads(cores)
* `-h`: print SAM header before reads
* `-f`: format of output file (default is SAM)
* `-F`: set [custom filter](https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax) - we will be using the filter to remove duplicates, multimappers and unmapped reads.

```bash
$ sambamba view -h -t 2 -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
wt_sample2_chip_sorted.bam > wt_sample2_chip_final.bam
```

We filter out unmapped reads by specifying in the filter `not unmapped`, and duplicates with `not duplicate`. Also, among the reads that are aligned, we filter out multimappers by specifying `[XS] == null`. 'XS' is a tag generated by Bowtie2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned and more than one alignment is found for the read.

Now that the alignment files contain only uniquely mapping reads, we are ready to perform peak calling!

> ### Filtering out Blacklisted Regions
> Although we do not perform this step, it is common practice to apply an additional level of filtering to our BAM files. That is, we remove alignments that occur with defined Blacklisted Regions. **We will filter out blacklist regions post-peak calling.**
> 
> Blacklisted regions represent artifact regions that tend to show artificially high signal (excessive unstructured anomalous reads mapping). These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats and typically appear uniquely mappable so simple mappability filters applied above do not remove them. The ENCODE and modENCODE consortia have compiled blacklists for various species and genome versions including human, mouse, worm and fly. These blacklisted regions (coordinate files) can be filtered out from our alignment files before proceeding to peak calling.
> 
> If we wanted to filter blacklist regions at this point in our workflow, we would use the following code:
> 
> ``` 
> # DO NOT RUN
> $ bedtools intersect -v -abam wt_sample2_chip_final.bam -b mm10-blacklist.v2.bed > wt_sample2_chip_final_blacklist_filtered.bam
> ```
> 
> _bedtools is a suite of tools that we will discuss in more detail in a later lesson when blacklist filtering is applied._

For CUT&RUN and ATAC-seq, there are additional parameters that you want to explore, and we list them below:

<details>
	<summary><b><i>How do the parameters change for CUT&RUN?</i></b></summary>
	<br>
	<p> Duplicate removal is an optional step in many CUT&RUN analysis approaches. The default is usually to keep duplicates, because CUT&RUN increases the likelihood of biological duplicates. More specifically, nuclease cleavage of chromatin by its stereotypical nature is influenced by conformation of chromatin and/or nuclease bias, and shorter DNA fragments also increase the likelihood of identical reads that are originated from different cells (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4). Therefore, we should remove duplicate with caution: assess the library complexity first, and then check if there are unreasonbly high amount of duplications. If not, and your experiment does not over-amplify, you might not want to remove the duplicate.
		
There are no mention of multi-mapper removal step in CUT&RUN analysis approaches neither. Maybe because Bowtie2 defaults to search for multiple alignments, yet only reports the best one (i.e. not using the `-k` option)
* Filtering of BAM files by fragment size.
   * After alignment, fragments can be divided into ≤ 120-bp and > 120-bp fractions. For transcription factors or proteins with an expected punctate binding profile, you can use the ≤ 120-bp fraction which is likely to contain binding sites. The range can be increased depending on the protein of interest, and alternatively BAM files without filtering can also be used. 
   * Done using sambamba: 

```bash
sambamba view --format \
  bam --nthreads 6 \
  -F "((template_length > 0 and template_length < 120) or (template_length < 0 and template_length > -120))" $file | samtools view -b > bams_sizeSelect/${s}-sizeSelect.bam
```
</p>
  
  For CUT&RUN, there are additional parameters that can be used. Here, we list options that have been reported by other groups. <b>It might not be neccessary to include any or all of these options.</b>  <i>We encourage you to explore the literature that resemble your research and method, and decide what is best for your data.</i>
		
* `--end-to-end`: An opposite option of `--local`. Bowtie2 will search for alignments involving all of the read characters. This is also called an "untrimmed" or "unclipped" alignment, and is only used when trimming is done prior to alignment.
* `--very-sensitive`: A preset option that results in slower running, but more sensitive and more accurate result.
* `--no-mixed`: Suppress unpaired alignments for paired reads. Otherwise, without this option, when Bowtie2 cannot find a concordant or discordant alignment for a pair, it tries to find alignments for the individual mates.
* `--no-discordant`: Suppress discordant alignments for paired reads. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints.
* `-I 10 -X 700`: For specifying the size range of inserts. In this example, 10-700 bp in length is used to ignore any remaining adapter sequence at the 3’ ends of reads during mapping.
* `--dovetail`: The term 'dovetailing' describes mates which extend past one another. It is unusual but is frequently encountered in CUT&RUN experiments. This option indicates that dovetailed alignments should be considered as concordant.</p>
	
</details>

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

