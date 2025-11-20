#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
/BiO/Share/Tools/bwa-0.7.18/bwa mem -M -t 10 -R '@RG\tID:12N\tPL:ILLUMINA\tLB:12N\tSM:12N\tCN:UNIST' -v 3 /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta /BiO/Teach/BME42401/Data/12N_R1.fastq.gz /BiO/Teach/BME42401/Data/12N_R2.fastq.gz | /BiO/Share/Tools/samtools-1.21/samtools view --bam --with-header --threads 10 --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --output ./12N.bam
