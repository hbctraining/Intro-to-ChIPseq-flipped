#! /bin/bash

for fq in ~/chipseq_workshop/raw_data/ko_sample1_chip.fastq.gz ~/chipseq_workshop/raw_data/ko_sample2_chip.fastq.gz
do

sbatch -p short -t 0-2:00 -n 6 --job-name chipseq-analysis -o %j.out -e %j.err \
--wrap="sh ~/chipseq_workshop/scripts/chipseq_analysis_on_input_file.sh $fq"

sleep 1	    # wait 1 second between each job submission
  
done
