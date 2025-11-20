#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610-01
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
/BiO/Share/Tools/bwa-0.7.18/bwa mem -M -t 10 -R '@RG\tID:12T\tPL:ILLUMINA\tLB:12T\tSM:12T\tCN:UNIST' -v 3 /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta /BiO/Teach/BME42401/Data/12T_R1.fastq.gz /BiO/Teach/BME42401/Data/12T_R2.fastq.gz | /BiO/Share/Tools/samtools-1.21/samtools view --bam --with-header --threads 10 --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --output ./12T.bam
/BiO/Share/Tools/bwa-0.7.18/bwa mem -M -t 10 -R '@RG\tID:TUMOR\tPL:ILLUMINA\tLB:TUMOR\tSM:TUMOR\tCN:UNIST' -v 3 /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta /BiO/Teach/BME42401/Data/TUMOR_R1.fastq.gz /BiO/Teach/BME42401/Data/TUMOR_R2.fastq.gz | /BiO/Share/Tools/samtools-1.21/samtools view --bam --with-header --threads 10 --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --output ./TUMOR.bam
