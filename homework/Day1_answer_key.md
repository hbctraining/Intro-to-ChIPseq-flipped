# Day 1 Answer Key

## QC using FASTQC

1. Once you have figured out what argument to use, run FastQC with 4 threads/cores.

```bash
$ fastqc -o ~/chipseq_workshop/results/fastqc/ -t 4 ~/chipseq_workshop/raw_data/wt*chip.fastq.gz
```

2. Do you notice any difference in running time? Does multi-threading speed up the run?

**Answer**: Yes, multi-threading speeds up the run. Also, the progress bar on the screen shows that two samples are processed simultaneously, which is another indication of multi-threading.



## Alignment using bowtie2

1. After your job has completed, check the resulting .out and .err files. What information do you obtain from each file?

**Answer**: The .out file is empty. The .err file contains the alignment statistics.

2. Take a quick peek at a sample BAM file using `samtools view`. Does the information you see line up with the fields we described above?

```bash
samtool view -h ~/chipseq_workshop/results/bowtie2/wt_sample2_chip.bam | less
```
**Answer**: Using the `less` command is helpful in allowing your to scroll through the file. You should see lines corresponding to each read and a short header section.

3. What is the alignment rate for the wt_sample2_chip? Do you think the alignment is good?

**Answer**: The alignment rate is 98.34%. The alignment is good.
