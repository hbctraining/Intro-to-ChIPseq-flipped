# Understanding chromatin biology using high throughput sequencing workshop

### Learning Objectives
* Understand the necessity for, and use of, the command line interface (bash) and HPC for analyzing high-throughput sequencing data.
* Understand best practices for designing a ChIP-seq / CUT&RUN / ATAC-seq experiment.
* Perform the steps involved in going from raw FASTQ files to peak calls for an individual sample.
* Review qualitative ways to assess peak calls and if they support the hypothesis

# Installations

***All:***

* [FileZilla Client](https://filezilla-project.org/download.php?type=client) (make sure you get ‘FileZilla Client')

***Mac users:***

* Plain text editor like [Sublime text](http://www.sublimetext.com/) or similar

***Windows users:***

* [GitBash](https://git-scm.com/download/win)
* Plain text editor like [Notepad++](http://notepad-plus-plus.org/) or similar

## Notes
* These materials focus on the use of local computational resources at Harvard, which are **only accessible to Harvard affiliates**
* Non-Harvard folks can [download the data]() and set up to work on their local clusters (with the help of local system administrators)

### Instructions for Harvard researchers with access to HMS-RC's O2 cluster

To run through the code in the lessons below, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).

1. Log in using `ssh ecommonsID@o2.hms.harvard.edu` and enter your password.
2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 1G /bin/bash` to get on a compute node _or as specified in the lesson_.
3. Proceed only once your command prompt has the word `compute` in it.
4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.
