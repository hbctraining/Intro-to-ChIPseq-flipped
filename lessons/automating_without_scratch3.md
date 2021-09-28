---
title: "Automating QC and Alignment"
author: "Meeta Mistry, Radhika Khetani, Jihe Liu, Mary Piper, Will Gammerdinger"
date: "Sept 15th, 2021"
---

Approximate time: 70 minutes

## Learning Objectives

* Create a shell script to streamline the analysis for QC and alignment 
* Utilize the SLURM job scheduler to run the workflow shell script on all samples
* Describe the difference between serial and parallel jobs

## Automating the ChIP-seq analysis path from sequence reads to BAM files

Once you have optimized tools and parameters using a single sample, you can write a script to run the whole (or a portion of the) workflow on multiple samples.

Writing a reusable shell script ensures that every sample is run with **the exact same parameters**, and helps to keep **track of all the tools and their versions**. The shell script is like a lab notebook; in the future, you (or your colleagues) can go back and check the workflow for methods and versions, making your work not only more efficient, but also reproducible.

Before we start the script, let's check how many cores our interactive session has by using `sacct`. 

```bash
$ sacct
```

We need to have an interactive session with 6 cores, if you already have one you are set. If you have a session with fewer cores then `exit` out of your current interactive session and start a new one with `-c 6`.

```bash
$ srun --pty -p short -t 0-12:00 -c 6 --mem 4G --reservation=HBC3 /bin/bash
```

### More Flexibility with variables

We can write a shell script with specific files as input, but to make it more flexible and efficient, it makes more sense to have the script accept those files as input when we run the script. To be able to provide one or more inputs to any shell script on the command prompt (instead of *hard coding* files in the script itself), we want to introduce you to the idea and implementation of **Positional Parameters**.

For example, we can refer to the components of the following command as numbered variables **within** the actual script:

```bash
# * DO NOT RUN *
$ sh run_analysis.sh input.fastq input.gtf 12
```

`$0` => run_analysis.sh

`$1` => input.fastq

`$2` => input.gtf

`$3` => 12

The variables $1, $2, $3,...$9 and so on are **positional parameters** in the context of the shell script, and can be used within the script to refer to the files/number specified on the command line. Basically, this example script is written with the expectation that $1 will be a fastq file and $2 will be a GTF file, and so on.

*There can be virtually unlimited numbers of inputs to a shell script, but it is wise to only have a few inputs to avoid errors and confusion when running a script that used positional parameters.*

> [This is an example of a simple script that used the concept of positional parameters and the associated variables](http://steve-parker.org/sh/eg/var3.sh.txt). You should try this script out after the class to get a better handle on positional parameters for shell scripting.

Let's use this new concept in the script we are creating. We want the first positional parameter ($1) to be the name of our fastq file. 

First, we need to start a new script called `chipseq_analysis_on_input_file.sh` in the `~/chipseq_workshop/scripts/` directory:

```bash
$ cd ~/chipseq_workshop/scripts/

$ vim chipseq_analysis_on_input_file.sh
```

When we run this script, we want to be able to say `sh chipseq_analysis_on_input_file.sh <name-of-fastq-file>`. To make this work, we have to replace all the places in the script where we want to refer to the fastq file, with the variable `$1`. However, this variable name is not intuitive, so we want to create a new variable called `fq` and copy the contents of `$1` into it.

```bash
#!/bin/bash/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
```

Now, when we need to use the fastq file in any command we can use `$fq` instead of the non-intuitive variable name `$1`.

> When we *set up variables* we do not use the `$` before it, but when we *use the variable*, we always have to have the `$` before it.
>
> For example: 
>
> initializing the `fq` variable => `fq=$1`
>
> using the `fq` variable => `fastqc $fq`

To ensure that all the output files from the workflow are properly named with sample IDs, we should extract the "base name" (or sample ID) from the name of the input file. There is a command called `basename` that will do that! 

We can save the output of the `basename` command and save it as a new variable called `$base` that we can use downstream.

```bash
# grab base of filename for naming outputs
base=`basename $fq .fastq.gz`          
```

> **What is `basename`?**
>
> The `basename` command: this command takes a path or a name, and trims away all the information **before** the last `\`. If you specify the string to clear away **at the end**, it will do that as well. 
> 
> In this case, if the variable `$fq` contains the path *"~/chipseq_workshop/raw_data/wt_sample2_chip.fastq.gz"*, `basename $fq .fastq.gz` will output "wt_sample2_chip".
>
> ***
>
> **What do the backticks used with `basename` (``) do?**
>   
> To assign the value of the `basename` command to the `base` variable, we encapsulate the `basename $fq .fastq.gz` command in backticks. This syntax is necessary for assigning the output of a given command to a given variable.

Next, we'll setup our directory structure for outputs using the `-p` option. This will make sure that `mkdir` will create the directory only if it does not already exist, and it will not throw an error if it does exist.

```
# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
mkdir -p ~/chipseq_workshop/results/fastqc
mkdir -p ~/chipseq_workshop/results/bowtie2
```



> **What do the `{}` used with `$base` in the file name variables do?**
>   
> To make sure that the output files are properly names with the sampleID contained in the `$base` variable, we use the variable in the new file name e.g. `~/chipseq/results/bowtie2/${base}_final.bam`. But, why does `$base` need to be inside the `{}`? 
>
> The curly brackets ensure that shell recognizes that only `base` is the variable name, and does not accidentally use the full string after the `$` as the variable name, e.g. `base_final.bam`. 

### Keeping track of tool versions

All of our variables are now staged. Next, let's make sure all the modules are loaded. This is also a good way to keep track of the versions of tools that you are using in the script:

```
# set up the software environment
module load fastqc/0.11.3
module load gcc/6.2.0
module load bowtie2/2.2.9
module load samtools/1.9
module load sambamba/0.7.1
```

### Preparing for future debugging

It is a good idea to use the `echo` command for debugging. `echo` prints the string of characters specified within the quotations onto the Terminal. When you strategically place `echo` commands indicating the step of the analysis, you can determine where the script failed using the last `echo` statement displayed when you are troubleshooting. For example, we could add the below `echo` command with `$base` before we perform the FastQC analysis to indicate the step and the sample ID.

```
echo "FastQC analysis of $base"
```

> You can also use `set -x`:
>
> `set -x` is a debugging tool that will make bash display the command before executing it. In case of an issue with the commands in the shell script, this type of debugging lets you quickly pinpoint the step that is throwing an error. Often, tools will display the error that caused the program to stop running, so keep this in mind for times when you are running into issues where this is not available.
> You can turn this functionality off by saying `set +x`

### Running the tools

Let's write up the commands to run the tools we have already tested, with a couple of modifications:
* use variable names instead of actual file names
* multithread when possible (`bowtie2`, `sambamba`)

```
# Run FastQC
fastqc -o $fastqc_out $fq

# Run bowtie2
bowtie2 -p 2 -q --local -x $genome -U $fq -S $align_sam

# Create BAM from SAM
samtools view -h -S -b -o $align_bam $align_sam

# Remove SAM file
rm $align_sam

# Sort BAM file by genomic coordinates
samtools sort $align_bam -o $align_sorted 

# Filter out multi-mappers and duplicates
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped " $align_sorted > $align_final

# Create indices for all the bam files for visualization and QC
samtools index $align_final
```

### Last addition to the script

It is best practice to have the script **usage** specified at the top of any script. This usage message should have information such that when your future self (or a co-worker) uses the script, they know what it will do and what input(s) are needed. For our script, we should have the following lines of comments right at the top after `#!/bin/bash/`:

```
# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Fastq files are aligned against the mm10 genome using Bowtie2. The outputted BAM file **does not** contain duplicate reads and is sorted by genomic coordinates using sambamba and samtools, respectively.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>
```

It is okay to specify this after all the commands are in the script, when you will have most clarity about the script's purpose.

Your script should now look like this:

```
#!/bin/bash/

# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Fastq files are aligned against the mm10 genome using Bowtie2. The outputted BAM file **does not** contain duplicate reads and is sorted by genomic coordinates using sambamba and samtools, respectively.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# grab base of filename for naming outputs
base=`basename $fq .fastq.gz`  

# directory with bowtie genome index
genome=/n/groups/shared_databases/bowtie2_indexes/mm10

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
mkdir -p ~/chipseq_workshop/results/fastqc
mkdir -p ~/chipseq_workshop/results/bowtie2

# set up output filenames and locations
fastqc_out=~/chipseq_workshop/results/fastqc
bowtie_results=~/chipseq_workshop/results/bowtie2

## set up file names
align_sam=~/chipseq_workshop/results/bowtie2/${base}.sam
align_bam=~/chipseq_workshop/results/bowtie2/${base}.bam
align_sorted=~/chipseq_workshop/results/bowtie2/${base}_sorted.bam
align_final=~/chipseq_workshop/results/bowtie2/${base}_final.bam

# set up the software environment
module load fastqc/0.11.3
module load gcc/6.2.0  
module load bowtie2/2.2.9
module load samtools/1.9
module load sambamba/0.7.1

echo "FastQC analysis of $base"

# Run FastQC and place the output to the appropriate folder
fastqc -o $fastqc_out $fq

echo "Bowtie2 alignment of $base"

# Run bowtie2
bowtie2 -p 2 -q --local -x $genome -U $fq -S $align_sam

echo "Convert SAM to BAM for $base"

# Create BAM from SAM
samtools view -h -S -b -o $align_bam $align_sam

# Remove SAM file
rm $align_sam

echo "Filtering $base BAM file"

# Sort BAM file by genomic coordinates
samtools sort $align_bam -o $align_sorted 

# Filter out duplicates
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped " $align_sorted > $align_final

# Create indices for all the bam files for visualization and QC
samtools index $align_final

echo "The analysis for $base is done!"
```

### Saving and running script

We should all have an interactive session with 6 cores, so we can run the script as follows:

```bash
$ sh chipseq_analysis_on_input_file.sh ~/chipseq_workshop/raw_data/wt_sample2_chip.fastq.gz
```

## Submitting jobs **in serial** to the SLURM scheduler

The above script will run in an interactive session **one file at a time**. But the whole point of writing this script was to run it on all files at once. How do you think we can do this?

To run the above script **"in serial"** for all of the files via the job scheduler, we can create a separate submission script that will need 2 components:

1. **SLURM directives** at the **beginning** of the script. This is to let scheduler know what resources we need in order to run our job on the compute node(s).
2. a **`for`** loop that iterates through and runs the above script for all the fastq files.

Below is what this second script (`chipseq_analysis_on_allfiles.slurm`) would look like **\[DO NOT RUN THIS\]**:

```
# **\[DO NOT RUN THIS\]**

#!/bin/bash

#SBATCH -p short 		# partition name
#SBATCH -t 0-2:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -n 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --job-name chipseq 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard error will be written

# this `for` loop, will take chip-seq fastq files as input and output filtered BAM files ready for peak calling.

for fq in ~/chipseq_workshop/raw_data/*.fastq.gz
do
  echo "running analysis on $fq"
  sh chipseq_analysis_on_input_file.sh $fq
done
```

**But, we don't want to run the analysis on these 8 samples one after the other,** we want to run them "in parallel" as 8 separate jobs! 

> **Note:** If you create and run the above script, or something similar to it, i.e. with SLURM directives at the top, you should give the script name `.run` or `.slurm` as the extension. This will make it obvious that it is meant to submit jobs to the SLURM scheduler. 
>
> To the run the above script, you would have used the following command: `sbatch chipseq_analysis_on_allfiles.slurm`.  

## Submitting jobs **in parallel** to the SLURM scheduler

Parallelization will save you a lot of time with real (large) datasets. To parallelize our analysis, we will still need to write a second script that will call the original script we just wrote. We will still use a `for` loop, but we will be creating a regular shell script and we will be specifying the SLURM directives a little differently as input to the `for` loop. 

Use `vim` to start a new shell script called `chipseq_run_allfiles.sh`: 

```bash
$ vim chipseq_run_allfiles.sh
```

Due to the space limit in our home directory, we will loop through only two `KO` samples. The command being submitted within the `for` loop is `sbatch` with SLURM directives specified on the same line (similar to what we have been doing for `srun`):

```bash
#! /bin/bash

for fq in ~/chipseq_workshop/raw_data/ko_sample1_chip.fastq.gz ~/chipseq_workshop/raw_data/ko_sample2_chip.fastq.gz
do

sbatch -p short -t 0-2:00 -n 6 --job-name chipseq-analysis -o %j.out -e %j.err \
--wrap="sh ~/chipseq_workshop/scripts/chipseq_analysis_on_input_file.sh $fq"

sleep 1	    # wait 1 second between each job submission
  
done
```
> Please note that after the `sbatch` directives, the command `sh ~/chipseq_workshop/scripts/chipseq_analysis_on_input_file.sh $fq` is in quotes.

Finally, you run the script as follows:
```bash
sh chipseq_run_allfiles.sh
```
What you should see on the output of your screen would be the 2 jobIDs that are returned from the scheduler for each of the jobs that your script submitted.

You can use `sacct <login_ID>` to check progress.

Don't forget about the `scancel` command, should something go wrong and you need to cancel your jobs.

> **NOTE:** All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons which are considered by the system administrators when picking one for a given HPC environment. 

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
