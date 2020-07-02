#!/bin/bash
# Author: Hailiang Mei
# Date: 2020/6/30
# E-mail: hailiang.mei@riken.jp

# general mapping for CUT&RUN
# Global variables
ref_genome="/home/janedoe/Mus_musculus/mm10/Bowtie2Index/genome"
adapter="AGATCGGAAGAGC"


# QC for raw data
for sample in seq_file_A seq_file_B
do
    trim_galore ${sample}.R1.fastq ${sample}.R2.fastq --quality 20 --phred33 -a ${adapter} --stringency 6 -e 0.1 --trim-n --basename ${sample} --paired --retain_unpaired & 
done

# bowtie2 mapping for CUT&RUN data
for sample in seq_file_A seq_file_B
do
    bowtie2 -q -N 1 -L 25 --maxins 800 --no-mixed --no-discordant -p 5 -x ${ref_genome} -1 ${sample}_val_1.fq -2 ${sample}_val_2.fq -S ${sample}.sam 2>log_${sample} &
done

# remove PCR dupliates and keep unique reads for CUT&RUN data
for sample in seq_file_A seq_file_B
do
    sambamba-0.7.1-linux-static markup -r -t 2 ${sample}.sorted.bam ${sample}.sorted.rmduplicate.bam
    sambamba-0.7.1-linux-static view -h -f sam -o ${sample}.sorted.rmduplicate.sam ${sample}.sorted.rmduplicate.bam
    grep -v 'XS:i' ${sample}.sorted.rmduplicate.sam > ${sample}.sorted.rmduplicate.unique.sam
    perl choose_singleton_and_filter.pl ${sample}.sorted.rmduplicate.unique.sam ${sample}.sorted.rmduplicate.unique.single.sam
    sambamba-0.7.1-linux-static view -h -S --format=bam -o ${sample}.sorted.rmduplicate.unique.single.bam ${sample}.sorted.rmduplicate.unique.single.sam
done

# merge replicates and generate bigwig file for CUT&RUN data
for sample in seq_file_A seq_file_B
do
    samtools merge -h ${sample}_replicate1.sorted.bam ${sample}_merge_replicates.bam ${sample}_replicate1.sorted.bam ${sample}_replicate2.sorted.bam
    sambamba-0.7.1-linux-static sort -o ${sample}_merge_replicates.sort.bam ${sample}_merge_replicates.bam
    sambamba-0.7.1-static index ${sample}_merge_replicates.sort.bam ${sample}_merge_replicates.sort.bam.bai
    bamCoverage --bam ${sample}_merge_replicates.sort.bam --outFileName ${sample}_merge_replicates.bigwig --outFileFormat bigwig --binSize 50 --numberOfProcessors 5 --normalizeUsing RPKM --extendReads 200 --ignoreDuplicates
done


# Allelic mapping for CUT&RUN
ref_genome="/home/janedoe/Mus_musculus/mm10/bowtie2_B6_PWK_PhJ_N_masked/genome"

# bowtie2 allelic mapping
for sample in seq_file_A seq_file_B
do
    bowtie2 -q -N 1 -L 25 --maxins 800 --no-mixed --no-discordant -p 5 -x ${ref_genome} -1 ${sample}_val_1.fq -2 ${sample}_val_2.fq 2>log_${sample} &
done

# snpsplit to determine the origin
for sample in seq_file_A
do
    SNPsplit --snp_file /home/janedoe/reference_genome/all_SNPs_PWK_GRCm38_merge.txt --paired --singletons -o SNPsplit_out/ --no_sort --conflicting ${sample}.sorted.bam &
done

# remove duplicates and keep unique
for sample in seq_file_A seq_file_B seq_file_C
do
    sambamba-0.7.1-linux-static markup -r -t 2 ${sample}.sorted.genome1.bam ${sample}.sorted.genome1.rmduplicate.bam
    sambamba-0.7.1-linux-static markup -r -t 2 ${sample}.sorted.genome2.bam ${sample}.sorted.genome2.rmduplicate.bam
    sambamba-0.7.1-linux-static view -h -f sam -o ${sample}.sorted.genome1.rmduplicate.sam ${sample}.sorted.genome1.rmduplicate.bam
    sambamba-0.7.1-linux-static view -h -f sam -o ${sample}.sorted.genome2.rmduplicate.sam ${sample}.sorted.genome2.rmduplicate.bam
    grep -v 'XS:i' ${sample}.sorted.genome1.rmduplicate.sam > ${sample}.sorted.genome1.rmduplicate.unique.sam
    grep -v 'XS:i' ${sample}.sorted.genome2.rmduplicate.sam > ${sample}.sorted.genome2.rmduplicate.unique.sam
    perl choose_singleton_and_filter.pl ${sample}.sorted.genome1.rmduplicate.unique.sam ${sample}.sorted.genome1.rmduplicate.unique.single.sam
    perl choose_singleton_and_filter.pl ${sample}.sorted.genome2.rmduplicate.unique.sam ${sample}.sorted.genome2.rmduplicate.unique.single.sam
    sambamba-0.7.1-linux-static view -h -S --format=bam -o ${sample}.sorted.genome1.rmduplicate.unique.single.bam ${sample}.sorted.genome1.rmduplicate.unique.single.sam
    sambamba-0.7.1-linux-static view -h -S --format=bam -o ${sample}.sorted.genome2.rmduplicate.unique.single.bam ${sample}.sorted.genome2.rmduplicate.unique.single.sam
    sambamba-0.7.1-linux-static markup -r -t 2 ${sample}.sorted.genome1_st.bam ${sample}.sorted.genome1_st.rmduplicate.bam
    sambamba-0.7.1-linux-static markup -r -t 2 ${sample}.sorted.genome2_st.bam ${sample}.sorted.genome2_st.rmduplicate.bam
    sambamba-0.7.1-linux-static view -h -f sam -o ${sample}.sorted.genome1_st.rmduplicate.sam ${sample}.sorted.genome1_st.rmduplicate.bam
    sambamba-0.7.1-linux-static view -h -f sam -o ${sample}.sorted.genome2_st.rmduplicate.sam ${sample}.sorted.genome2_st.rmduplicate.bam
    grep -v 'XS:i' ${sample}.sorted.genome1_st.rmduplicate.sam > ${sample}.sorted.genome1_st.rmduplicate.unique.sam
    grep -v 'XS:i' ${sample}.sorted.genome2_st.rmduplicate.sam > ${sample}.sorted.genome2_st.rmduplicate.unique.sam
    sambamba-0.7.1-linux-static view -h -S --format=bam -o ${sample}.sorted.genome1_st.rmduplicate.unique.bam ${sample}.sorted.genome1_st.rmduplicate.unique.sam
    sambamba-0.7.1-linux-static view -h -S --format=bam -o ${sample}.sorted.genome2_st.rmduplicate.unique.bam ${sample}.sorted.genome2_st.rmduplicate.unique.sam
done

# merge and generate bigwig files
for sample in seq_file_A seq_file_B
do
    samtools merge -h ${sample}.sorted.genome1.rmduplicate.unique.single.bam ${sample}.genome1.merge.allelic.bam ${sample}.sorted.genome1.rmduplicate.unique.single.bam ${sample}.sorted.genome1_st.rmduplicate.unique.bam
    samtools merge -h ${sample}.sorted.genome2.rmduplicate.unique.single.bam ${sample}.genome2.merge.allelic.bam ${sample}.sorted.genome2.rmduplicate.unique.single.bam ${sample}.sorted.genome2_st.rmduplicate.unique.bam
    sambamba-0.7.1-linux-static sort -o ${sample}.genome1.merge.allelic.sort.bam ${sample}.genome1.merge.allelic.bam
    sambamba-0.7.1-linux-static sort -o ${sample}.genome2.merge.allelic.sort.bam ${sample}.genome2.merge.allelic.bam
    bamCoverage --bam ${sample}.genome1.merge.allelic.sort.bam --outFileName ${sample}.genome1.merge.allelic.sort.bigwig --outFileFormat bigwig --binSize 50 --numberOfProcessors 5 --normalizeUsing RPKM --extendReads 200 --ignoreDuplicates
    bamCoverage --bam ${sample}.genome2.merge.allelic.sort.bam --outFileName ${sample}.genome2.merge.allelic.sort.bigwig --outFileFormat bigwig --binSize 50 --numberOfProcessors 5 --normalizeUsing RPKM --extendReads 200 --ignoreDuplicates

done
