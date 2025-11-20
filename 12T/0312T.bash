#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610-03
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
/BiO/Share/Tools/gatk-4.6.1.0/gatk Mutect2 --java-options "-XX:+UseSerialGC -Xmx20g" --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --input 12T.Sort.MarkDuplicates.BQSR.bam --input TUMOR.Sort.MarkDuplicates.BQSR.bam --12T-sample 12T --output SAMPLE.vcf --native-pair-hmm-threads 10 --max-mnp-distance 0
/BiO/Share/Tools/gatk-4.6.1.0/gatk FilterMutectCalls --java-options "-XX:+UseSerialGC -Xmx20g" --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --variant SAMPLE.vcf --output SAMPLE.filter.vcf
/usr/bin/awk -F '	' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' SAMPLE.filter.vcf > SAMPLE.PASS.vcf
/BiO/Share/Tools/gatk-4.6.1.0/gatk IndexFeatureFile --java-options "-XX:+UseSerialGC -Xmx20g" --input SAMPLE.PASS.vcf --output SAMPLE.PASS.vcf.idx
/usr/bin/perl /BiO/Share/Tools/vcf2maf-1.6.21/vcf2maf.pl --vep-path /BiO/Share/Tools/ensembl-vep-release-110.1/ --vep-data /BiO/Share/Tools/ensembl-vep-release-110.1/ --vep-forks 10 --ncbi-build 'GRCh38' --input-vcf SAMPLE.PASS.vcf --output SAMPLE.PASS.maf --tumor-id TUMOR --12T-id 12T --ref-fasta /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --vep-overwrite
