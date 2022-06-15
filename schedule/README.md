# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to the command-line interface](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/) workshop is a prerequisite.


## Pre-reading:

* Please **study the contents** and **work through all the exercises** within the following lessons:
  * [Shell basics review](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/shell_review.html)
  * [Best Practices in Research Data Management (RDM)](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/04a_data_organization.html)
  
  
## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop Introduction](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lectures/Intro_to_workshop.pdf) | Meeta |
| 09:45 - 11:00 | [Introduction to ChIP-seq](https://www.dropbox.com/s/i7m8a95sasoerwv/Introduction%20to%20ChIP-seq%202021.pdf?dl=1) | Meeta |
| 11:00- 11:05 | Break|  |
| 11:05 - 11:50 | [Working in an HPC environment](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/03_working_on_HPC.html) | Jihe |
| 11:50 - 12:00 | Overview of self-learning materials and homework submission | Jihe |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Experimental design considerations and understanding the ChIP-seq workflow](../lessons/01_ChIPseq_design_and_workflow.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Before you begin thinking about performing a ChIP-seq experiment, it is important to plan for it. There are many things to consider depending on the cells you are working with, and your protein of interest. <br><br>In this lesson, we will:<br>
             - Describe different types of controls and how they can help<br>
             - Highlight sequencing considerations different binding profiles<br>
             - Introduce you to the ChIP-seq analysis workflow R<br><br>
        </details>
   
   2. [Dataset overview and project organization](../lessons/02_dataset_and_project_setup.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>We are ready (and excited!) to get started with our ChIP-seq analysis. But, wait! There are just a few things to do before we get our hands on the data. <br><br>In this lesson you will:<br>
            - Learn about the dataset we are using in this workshop<br>
            - Organize your space, so you are setup for success<br><br>
         </details>

   3. [Quality Control of Sequence Data: Running FASTQC and evaluating results](../lessons/03_QC_FASTQC.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>The first step of most NGS analyses is to evaluate the quality of your sequencing reads. <br><br>In this lesson you will explore:<br>
            - The FASTQC software, and how to run it on your raw sequencing data<br>
            - The HTML report that is returned from FASTQC and how to interepret the different plotsn<br><br>
        </details>
        
   4. [Alignment using Bowtie2](../lessons/04_alignment_using_bowtie2.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>.  ...... <br><br>In this lesson you will learn:<br>
            - Reading different types (formats) of data<br>
            - Inspecting the contents and structure of the dataset once you have read it in<br><br>
        </details>


    > **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).
    > 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (replace the "XX" in the username with the number you were assigned in class). 
    > 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 1G /bin/bash` to get on a compute node or as specified in the lesson.
    > 3. Proceed only once your command prompt has the word `compute` in it.
    > 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox](https://www.dropbox.com/request/czRTtxT1fAHN1Ya6mIw7) the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:15 | Self-learning lessons review |  All |
| 10:15 - 11:00 | [Filtering BAM files](../lessons/05_filtering_BAM_files.md) | Jihe |
| 11:00- 11:05 | Break|  |
| 11:05 - 12:00 | [Peak calling](../lessons/06_peak_calling_macs.md) | Meeta |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:

 * [Handling peak files using `bedtools`](../lessons/07_handling_peaks_bedtools.md)
 * [Creating bigWig files](../lessons/08_creating_bigwig_files.md)
 * [Qualitative assessment of peak enrichment using deepTools](../lessons/09_data_visualization.md)

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox](https://www.dropbox.com/request/HOtwHGI6pH9Ha8txuvo3) the **day before the next class**.
   
### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:00 | Self-learning lessons review | All |
| 10:00 - 10:30 | [Troubleshooting your ChIP-seq analysis](../lessons/troubleshooting_chipseq_partI.md) | Meeta |
| 10:30- 10:35 | Break|  |
| 10:35 - 11:45 | [Automating the ChIP-seq workflow](../lessons/10_automation_new.md) | Jihe |
| 11:45 - 12:00 | [Wrap-up](../lectures/Wrap-up_new.pdf) | Meeta |


## Answer keys

* Day 1 exercises 
  * [Data Management and project organization](../homework/Day1_readme_answerkey.md)
  * [QC and Alignment questions](../homework/Day1_answer_key.md)

* Day 2 exercises 
  * [Handling peak calls](../homework/Day2_answer_key.md)

* Day 3 In-class 
  * [Automation shell script](../homework/chipseq_analysis_on_input_file.sh)
  * [Parallelization script](../homework/chipseq_run_allfiles.sh)

***


## Building on this workshop
* [Integration of ChIP-seq and RNA-seq](../lessons/integrating_rna-seq_and_chip-seq.md)
* [Advanced bash commands (aliases, copying files, and symlinks)](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/more_bash_cluster.html)
* [Introduction to R workshop materials](https://hbctraining.github.io/Intro-to-R-flipped/#lessons) 


***
*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
