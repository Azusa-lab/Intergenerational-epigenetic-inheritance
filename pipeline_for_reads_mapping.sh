#!/bin/bash
# Author: Hailiang Mei
# Date: 2020/6/30
# E-mail: hailiang.mei@riken.jp

# Global variables
ref_genome="/home/janedoe/Mus_musculus/mm10/Bowtie2Index/genome"
adapter="AGATCGGAAGAGC"


# QC for raw data
for sample in seq_file_A seq_file_B
do
    trim_galore ${sample}.R1.fastq ${sample}.R2.fastq --quality 20 --phred33 -a ${adapter} --stringency 6 -e 0.1 --trim-n --basename ${sample} --paired --retain_unpaired & 
done

# bowtie2 mapping
for sample in seq_file_A seq_file_B
do
    bowtie2 -q -N 1 -L 25 --maxins 800 --no-mixed --no-discordant -p 5 -x ${ref_genome} -1 ${sample}_val_1.fq -2 ${sample}_val_2.fq -S ${sample}.sam 2>log_${sample} &
done

# remove PCR dupliates
