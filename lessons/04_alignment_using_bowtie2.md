---
title: "Alignment using Bowtie2"
author: "Mary Piper, Radhika Khetani, Jihe Liu"
date: "Aug 17th, 2021"
---

Contributors: Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu

Approximate time:

**Link to issue describing the modifications to be made:** https://github.com/hbctraining/Intro-to-ChIPseq-flipped/issues/7

## Learning Objectives

* Learn how to align reads to the genome using Bowtie2
* Understand SAM and BAM file format, and learn how to convert SAM file to BAM file
* Run alignment script and evaluate the alignment result

## Alignment to Genome

Now that we have assessed the quality of our sequence data, we are ready to align the reads to the reference genome. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) is a fast and accurate alignment tool that indexes the genome with an FM Index based on the Burrows-Wheeler Transform method to keep memory requirements low for the alignment process. *Bowtie2* supports gapped, local and paired-end alignment modes, and it works best for reads that are at least 50 bp (shorter read lengths should use Bowtie1). By default, Bowtie2 performs a global end-to-end read alignment, which is best for quality-trimmed reads. However, it also offers a local alignment mode, which will perform soft-clipping for the removal of poor quality bases or adapters from untrimmed reads. We will use this option since we did not trim our reads.

<p align="center">
<img src="../img/chip_workflow_june2017_step1_align.png" width="400">
</p>

> #### How do other aligners compare?
> We use Bowtie2 to align our reads in this workshop, but there are a number of other options. For **[bwa](http://bio-bwa.sourceforge.net/)**, the mapping rates are higher, with an equally similar increase in the number of duplicate mappings. Consequently, there is a significantly higher number of mapped reads and a much larger number of peaks being called (30% increase compared to Bowtie2). When we compare the peak calls generated from different aligners, the **bwa** peak calls are a superset of those called from the Bowtie2 aligments. It is yet to be determined whether or not these additional peaks are true positives. 

### Bowtie2 index

To perform the Bowtie2 alignment, a genome index is required. The index is analagous to the index in a book. By indexing the genome, we have organized it in a manner that now allows for efficient search and retrieval of matches of the query (sequence read) to the genome.

The O2 cluster has a designated directory (`/n/groups/shared_databases/`) with shared databases for human and other commonly used model organisms, where O2 users could readily access. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool-specific data, and data from public databases, such as NCBI and PDB. Therefore, when you use a tool that requires a reference, it is worth taking a quick look at this directory and checking whether your desired reference is already deposited here.

>```bash
>$ ls -l /n/groups/shared_databases/
>```

Back to our data, we will use mouse `mm10` version as the reference genome. The index is located at `/n/groups/shared_databases/bowtie2_indexes/mm10`, so we will directly refer to it. However, if you need to create a genome index yourself, you could use the following command:

```bash
# DO NOT RUN

bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>
```

### Aligning reads to the genome with Bowtie2

We could now start with the read alignment. Let's first create a `bowtie2` directory:

```bash
# Create bowtie2 directory
mkdir ~/chipseq_workshop/results/bowtie2
```

We then need to load the module. We could find out more about bowtie2 on O2:

```bash
$ module spider bowtie2
```

Notice that before we load bowtie2, we also need to load the gcc compiler (as is the case for many other NGS analysis tools on O2). As a tip, we recommend always run `module spider` first to check any dependent modules.

```bash
$ module load gcc/6.2.0 bowtie2/2.2.9
```

The below is an example code to run bowtie2 on a single FASTQ file `wt_sample2_chip`. Details on Bowtie2 and its functionality can be found in the [user manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml); we encourage you to peruse through to get familiar with all available options.

```bash
# DO NOT RUN
$ bowtie2 -p 2 -q --local \
-x /n/groups/shared_databases/bowtie2_indexes/mm10 \
-U ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz \
-S ~/chipseq_workshop/results/wt_sample2_chip.sam
```

Some basic options for aligning reads to the genome using Bowtie2 are:

* `-p`: number of processors/cores
* `-q`: reads are in FASTQ format
* `--local`: local alignment feature to perform soft-clipping
* `-x`: /path/to/genome_indices_directory
* `-U`: /path/to/FASTQ_file
* `-S`: /path/to/output/SAM_file

## Alignment file format: SAM/BAM

The output from the Bowtie2 aligner is an unsorted SAM file, also known as **Sequence Alignment Map format**. The SAM file is a **tab-delimited text file** that contains information for each individual read and its alignment to the genome. While we will go into some features of the SAM format, the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Each section begins with character ‘@’ followed by **a two-letter record type code**. These are followed by two-letter tags and values. Example of some common sections are provided below:

```
@HD  The header line
VN: format version
SO: Sorting order of alignments

@SQ  Reference sequence dictionary
SN: reference sequence name
LN: reference sequence length
SP: species

@PG  Program
PN: program name
VN: program version
```

Following the header is the **alignment section**. Each line corresponds to the alignment information for a single read. Each alignment line has **11 mandatory fields for essential mapping information** and a variable number of other fields for aligner-specific information. 

![SAM1](../img/sam_bam.png)

An example read mapping is displayed above. *Note that the example above spans two lines, but in the actual file it is a single line.* Let's go through the fields one at a time. 

- **`QNAME`:** Query name or read name - this is the same read name present in the header of the FASTQ file
- **`FLAG`:** numerical value providing information about read mapping and whether the read is part of a pair.
 
  > **NOTE:** The information stored inside the FLAG is additive based on the following information being TRUE or FALSE:
  > 
  > | Flag | Description |
  > | ------:|:----------------------:|
  > | 1 | read is mapped |
  > | 2 | read is mapped as part of a pair |
  > | 4 | read is unmapped |
  > | 8 | mate is unmapped |
  > | 16| read reverse strand|
  > | 32 | mate reverse strand |
  > | 64 | first in pair |
  > | 128 | second in pair |
  > | 256 | not primary alignment |
  > | 512 | read fails platform/vendor quality checks |
  > | 1024| read is PCR or optical duplicate |
  > 
  > * For a given alignment, each of these flags are either **on or off**, indicating the condition is **true or false**. 
  > * The `FLAG` is a combination of all of the individual flags (from the table above) that are true for the alignment 
  > * The beauty of the flag values is that **any combination of flags can only result in one sum**.
  > 
  > **There are tools that help you translate the bitwise flag, for example [this one from Picard](https://broadinstitute.github.io/picard/explain-flags.html)**

- **`RNAME`:** is the reference sequence name, giving the chromosome to which the read maps. The example read is from chromosome 1, which explains why we see 'chr1'. 
- **`POS`:** refers to the 1-based leftmost position of the alignment. 
- **`MAPQ`:** is giving us the alignment quality, the scale of which will depend on the aligner being used. 
- **`CIGAR`:** is a sequence of letters and numbers that represent the *edits or operations* required to match the read to the reference. The letters are operations that are used to indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each 'operation'.

Now to the remaining fields in our SAM file:

![SAM1](../img/sam_bam3.png)

The next three fields are more pertinent to paired-end data. 

- **`MRNM`:** is the mate reference name. 
- **`MPOS`:** is the mate position (1-based, leftmost). 
- **`ISIZE`:** is the inferred insert size.

Finally, you have the raw sequence data from the original FASTQ file stored for each read:

- **`SEQ`:** is the raw sequence
- **`QUAL`:** is the associated quality values for each position in the read.

Since we haven't run the code yet, we could take a quick peek at a sample SAM file here at `/n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip.sam`. Since it is just a text file, we can browse through it using `less`:

``` bash
$ less /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip.sam
```

**Does the information you see line up with the fields we described above?**

### Changing file format from SAM to BAM

While the SAM alignment file from Bowtie2 is human readable, we need a BAM alignment file for downstream analysis. A BAM file is a binary equavalent version of SAM file, in other words, the same file in a compressed format. Therefore, BAM file is not human readable, and it is much smaller in size. BAM file is the typical format used in bioinformatics tools. We will use [Samtools](http://samtools.github.io) to convert the file format from SAM to BAM.

> NOTE: Once we generate the BAM file, we don't need to retain the SAM file anymore - we could delete it to save space.

Let's load the module `samtools`:

```bash
$ module load gcc/6.2.0 # you may not need to load this if you are working in the same session from Bowtie2
$ module load samtools/1.9
```

You can find detailed instructions for different samtools functions in this [manual](http://www.htslib.org/doc/samtools-1.2.html). For our purpose, we will use the command `samtools view` with the following parameters:

* `-h`: include header in output
* `-S`: input is in SAM format
* `-b`: output BAM format
* `-o`: /path/to/output/file

```bash
# DO NOT RUN
$ samtools view -h -S -b \
-o ~/chipseq_workshop/results/wt_sample2_chip.bam \
~/chipseq_workshop/results/wt_sample2_chip.sam
```

## Create SBATCH script for the alignment

Genome alignment usually takes quite a while to finish - that's why we don't run the codes on an interactive node. Instead, we will create a SBATCH script, `alignment.sbatch` under the `~/chipseq_workshop/` directory, and submit this script as a job on the cluster. Let's specify the job submission options as below (don't forget the shebang line, `#!/bin/bash` at the begining):

```
#SBATCH -p short              # partition name
#SBATCH -c 2                  # number of cores
#SBATCH -t 0-2:00             # time limit
#SBATCH --mem 8G              # requested memory
#SBATCH --job-name alignment 	# job name
#SBATCH -o %j.out			          # file to which standard output will be written
#SBATCH -e %j.err 		          # file to which standard error will be written
```

In the body of the script, we will load the required modules, run bowtie2 to obtain alignment SAM file, and then convert SAM file to BAM file using samtools. Please refer to the corresponding codes we discussed earlier in this lesson, to come up with the whole script. Once you are done, submit the script as a job, using `sbatch alignment.sbatch` command.

<details>
  <summary>Solution</summary>
 
``` bash
#!/bin/bash
 
#SBATCH -p short              # partition name
#SBATCH -c 2                  # number of cores
#SBATCH -t 0-2:00             # time limit
#SBATCH --mem 8G              # requested memory
#SBATCH --job-name alignment 	# job name
#SBATCH -o %j.out			          # file to which standard output will be written
#SBATCH -e %j.err 		          # file to which standard error will be written

module load gcc/6.2.0 bowtie2/2.2.9 samtools/1.9
 
bowtie2 -p 2 -q --local \
-x /n/groups/shared_databases/bowtie2_indexes/mm10 \
-U ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz \
-S ~/chipseq_workshop/results/wt_sample2_chip.sam
 
samtools view -h -S -b \
-o ~/chipseq_workshop/results/wt_sample2_chip.bam \
~/chipseq_workshop/results/wt_sample2_chip.sam

rm ~/chipseq_workshop/results/wt_sample2_chip.sam        
```

</details>

> NOTE:
> - The job takes about 50 minutes to finish. You could monitor the progress using the `sacct` command;
> - In the last line of the solution code, we remove the SAM file after generating the BAM file. We recommend do so to save space.

**Exercise:**

After the running is finished, check the resulting `.out` and `.err` files. What information do you obtain from each file? WHat is the alignment rate for the `wt_sample2_chip`? Do you think the alignment is good?

