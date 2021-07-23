# lesson 2: QC
- Code
```
# Log into compute node
srun --pty -p interactive -t 0-12:00 --mem 8G /bin/bash

# Create the workshop project
cp -R /n/groups/hbctraining/harwell-datasets/chipseq_workshop ~

# Load module
module load fastqc/0.11.3

# Run fastqc. Take 2 min
fastqc -o ~/chipseq_workshop/results/fastqc/ ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz

# Log into compute node with multiple cores
srun --pty -c 4 -p interactive -t 0-12:00 --mem 8G /bin/bash

# Run fastqc again. Take less time
fastqc -o ~/chipseq_workshop/results/fastqc/ -t 4 ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz
```

- Note
1. The qualities of all samples (WT and KO) are good.
2. Review in detail one of the FASTQC reports

# lesson 3: Alignment;
- Code
```
# Create bowtie2 directory
mkdir ~/chipseq_workshop/results/bowtie2
```

> example code to show in the material
```
bowtie2 -p 2 -q --local \
-x /n/groups/shared_databases/bowtie2_indexes/mm10 \
-U ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz \
-S ~/chipseq_workshop/results/wt_sample2_chip.sam

samtools view -h -S -b \
-o ~/chipseq_workshop/results/wt_sample2_chip.bam \
~/chipseq_workshop/results/wt_sample2_chip.sam
```

> Exercise: run the alignment using sbatch file (the expected run time is 48 min)
```
#!/bin/sh

#SBATCH -p priority
#SBATCH -c 2
#SBATCH -t 0-2:00
#SBATCH --mem 8G

bowtie2 -p 2 -q --local \
-x /n/groups/shared_databases/bowtie2_indexes/mm10 \
-U ~/chipseq_workshop/data/wt_sample2_chip.fastq.gz \
-S ~/chipseq_workshop/results/wt_sample2_chip.sam

samtools view -h -S -b \
-o ~/chipseq_workshop/results/wt_sample2_chip.bam \
~/chipseq_workshop/results/wt_sample2_chip.sam
```

- Note
1. The rep1_chip WT sample has low alignment rate (35%). Other samples have > 98% alignment.
2. Present the code for bowtie2 (DO NOT RUN), discuss the parameters
3. Present code for samtools; discuss samtools and also the need for smaller file 
4. This lesson will include an sbatch script (the first one) as an exercise. Run one sample only as is currently in the lesson
> Another exercise question - look at your .err file - comment on the alignment rate
> Have a note that we are not actually using the SAM file, and we could delete it

# lesson 4: filtering (in class)
- Code
```
sambamba sort -t 2 \
    -o ~/chipseq_workshop/results/wt_sample2_chip_sorted.bam \
    ~/chipseq_workshop/results/wt_sample2_chip.bam

sambamba view -h -t 2 -f bam \
    -F "[XS] == null and not unmapped  and not duplicate" \
    ~/chipseq_workshop/results/wt_sample2_chip_sorted.bam > ~/chipseq_workshop/results/wt_sample2_chip_final.bam

# After the alignment, also need to index the final bam files.
samtools index ~/chipseq_workshop/results/wt_sample2_chip_final.bam
```

- Note
1. Running time: 7 min
 
# lesson 5: peak calling (in class)

- Code

> We will have them grab the snapshot files to run this since we have not created filtered BAM files for all 

```
# Create macs2 directory in results/
mkdir ~/chipseq_workshop/results/macs2

macs2 callpeak -t /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_chip_final.bam \
    -c /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_input_final.bam \
    -f BAM -g mm \
    -n wt_sample1 \
    --outdir ~/chipseq_workshop/results/macs2 2> ~/chipseq_workshop/results/macs2/wt_sample1_macs2.log

macs2 callpeak -t /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip_final.bam \
    -c /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_input_final.bam \
    -f BAM -g mm \
    -n wt_sample2 \
    --outdir ~/chipseq_workshop/results/macs2 2> ~/chipseq_workshop/results/macs2/wt_sample2_macs2.log
```

- Note:
1. Running time: 5 min (sample1), 6 min (sample2)
> In class exercise (?) Have them run this on the KO samples too? (but the results are not useful for later wt vs ko plotting, because only peak file for wt will be used).

# lesson 6: handling replicate (bedtools)

- Code
1. Filter out black listed region

```
bedtools intersect \
-v \ 
-a ~/chipseq_workshop/results/macs2/wt_sample1_peaks.narrowPeak \
-b ~/chipseq_workshop/reference/mm10-blacklist.v2.bed \
> ~/chipseq_workshop/results/macs2/wt_sample1_peaks_filtered.bed

bedtools intersect \
-v \ 
-a ~/chipseq_workshop/results/macs2/wt_sample2_peaks.narrowPeak \
-b ~/chipseq_workshop/reference/mm10-blacklist.v2.bed \
> ~/chipseq_workshop/results/macs2/wt_sample2_peaks_filtered.bed
```

2. Finalize overlap peaks with more stringent criteria

```
# Final code to generate the result
bedtools intersect -a ~/chipseq_workshop/results/macs2/wt_sample1_peaks_filtered.bed -b ~/chipseq_workshop/results/macs2/wt_sample2_peaks_filtered.bed -wo -f 0.3 -r > ~/chipseq_workshop/results/macs2/wt_peaks_final.bed
```

- Note: 
1. Number of peaks before and after filtering out black listed region
WT_pair2: 14294 -> 14134
WT_pair3: 20013 -> 19829

# lesson7: create bigwigs
- Code
```
# Create visualization/bigWig directory in results directory
mkdir ~/chipseq_workshop/results/visualization/bigWig

# Create bigwig file
bamCoverage -b ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20
```

- Note
1. Create bigWig for one sample (wt sample2). And then discuss other options for bigWig including normalization, bamCompare.

# lesson8: Data visualization

> * Interactively compute matrix for both WT reps. 
> * For the ENCODE, provide code and explain it - but we will provide snapshot of compute matrix files
> * Total of 4 plots in this lesson: Wt reps, Encode Fig6A, Encode Fig6B, Wt versus KO
> * Exercise have them run compute Matrix and plotProfile for Wtrep2 and KOrep?

2. Create signal profile
- Plot1: WT reps
Compute matrix (**run 12 min in an interactive node**)

```
computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S ~/chipseq_workshop/results/visualization/bigWig/wt_sample1_chip.bw ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-p 6
```

Plot density
```
# Create figures directory under visualization
mkdir ~/chipseq_workshop/results/visualization/figures

plotProfile -m ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot1_wt.png \
--regionsLabel "" \
--perGroup \
--colors red blue \
--samplesLabel "PRDM16_sample1" "PRDM16_sample2" \
--refPointLabel "PRDM16 binding sites"
```

- Plot2: encode fig 6a
> NOTE: the prompt will show the message "The following chromosome names did not match between the bigwig files", which is fine. It is not an error message.
> **run 20 min in an interactive node**
```
computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw /n/groups/hbctraining/harwell-datasets/encode-chipseq/H3k04me1UE14_mm10.bw /n/groups/hbctraining/harwell-datasets/encode-chipseq/H3k04me3UE14_mm10.bw /n/groups/hbctraining/harwell-datasets/encode-chipseq/H3k27me3UE14_mm10.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_encode_matrix.gz \
-p 6
```

Plot density
```
plotProfile -m ~/chipseq_workshop/results/visualization/wt_encode_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot2_wt_encode.png \
--regionsLabel "" \
--perGroup \
--colors blue green red orange \
--samplesLabel "PRDM16" "H3K4me" "H3K4me3" "H3K27me3" \
--refPointLabel "PRDM16 binding sites"
```

- Plot3: encode fig 6b
```
computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw /n/groups/hbctraining/harwell-datasets/encode-chipseq/H3k27acUE14_mm10.bw /n/groups/hbctraining/harwell-datasets/encode-chipseq/CbellumH3k27acMAdult8wks_mm10.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_encode_acetylation_matrix.gz \
-p 6
```

Plot density
```
plotProfile -m ~/chipseq_workshop/results/visualization/wt_encode_acetylation_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot3_wt_encode_acetylation.png \
--regionsLabel "" \
--perGroup \
--colors blue red green \
--samplesLabel "PRDM16" "H3K27ac Embryonic" "H3K27ac Adult" \
--refPointLabel "PRDM16 binding sites"
```

- Plot4: WT vs KO
Script to process KO samples: provide snapshot of KO bigwig file
```
# Navigate to results directory
computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw /n/groups/hbctraining/harwell-datasets/workshop_material/results/visualization/bigWig/ko_sample2_chip.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_ko_matrix.gz \
-p 6

plotProfile -m ~/chipseq_workshop/results/visualization/wt_ko_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot4_wt_ko.png \
--regionsLabel "" \
--perGroup \
--colors blue red \
--samplesLabel "PRDM16_WT" "PRDM16_KO" \
--refPointLabel "PRDM16 binding sites"
```

- Note
1. Use the range of 4 kb before and after the center
2. Need to lift the encode profile from mm9 to mm10 (using CrossMap)

