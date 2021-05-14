# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to the command-line interface](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/) workshop is a prerequisite.


## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 10:00 - 10:15 | [Workshop Introduction]() |  |
| 10:15 - 11:15 | [Introduction to ChIP-seq]() | [Dr. Shannan Ho Sui]() |
| 11:15- 11:20 | Break|  |
| 11:20 - 11:45 | [Understanding the ChIP-seq workflow]() | |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | |

### Before the next class:

* Please **study the contents** and **work through all the exercises** within the following lessons:
  * [Shell basics review](../lessons/shell_review.md)
  * [Best Practices in Research Data Management (RDM)]()

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:00 - 9:30 | [Working in an HPC environment - Review]() |  |
| 10:20 - 11:05 | [Project Organization (using Data Management best practices)]() |  |
| 11:05- 11:10 | Break|  |
| 11:10 - 11:45 | [Quality Control of Sequence Data: Running FASTQC]() |  |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission |  |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:

 * [Sequence alignment and considerations for ChIP-seq data]()
 * [Aligning reads using Bowtie2]()
 * [Filtering BAM files]()

    > **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).
    > 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (replace the "XX" in the username with the number you were assigned in class). 
    > 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 1G /bin/bash` to get on a compute node or as specified in the lesson.
    > 3. Proceed only once your command prompt has the word `compute` in it.
    > 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox]() the **day before the next class**.
   
### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning lessons review | All |
| 10:30- 10:35 | Break|  |
| 10:35 - 11:35 | [Peak Calling with MACS2]() | |
| 11:35 - 12:00 | Overview of self-learning materials and homework submission |  |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:

 * [Creating bigWig files)]()
 * [Visualization of signal profiles]()

     > **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).
     > 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (replace the "XX" in the username with the number you were assigned in class). 
     > 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 8G /bin/bash` to get on a compute node or as specified in the lesson.
     > 3. Proceed only once your command prompt has the word `compute` in it.
     > 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox](https://www.dropbox.com/request/4a1W0PltytfzKKhl2Xax) the **day before the next class**.
   
### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 4

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:10 | Self-learning lessons review | All |
| 10:10 - 10:45 | [Troubleshooting ChIP-seq Data Analysis]() |  |
| 10:45 - 11:45 | [Automating the ChIP-seq workflow]()|  |
| 11:45 - 12:00 | [Wrap up](f) |  |

***

* Downloadable Answer Keys (Day 2 exercises): 
  * ...

* Downloadable Answer Keys (Day 3 exercises): 
  * ...

* [Automation Script])

***

## Resources

***

## Building on this workshop

***
*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
