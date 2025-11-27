#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=week08
#SBATCH --mem=20G
#SBATCH --export=ALL
cnvkit.py batch --processes 10 --normal "$(realpath ../week02/PS-826-B_DNA.Sort.MarkDuplicates.BQSR.bam)" --fasta "/BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta" --targets "/BiO/Share/Database/UCSC/hg38/UCSC_hg38.bed" --annotate "/BiO/Share/Database/UCSC/hg38/refFlat.txt"  --scatter --diagram --output-dir "$(realpath .)" "$(realpath ../week02/PS-826-T_DNA.Sort.MarkDuplicates.BQSR.bam)"
grep --invert-match '_' "$(realpath ./PS-826-T_DNA.Sort.MarkDuplicates.BQSR.cnr)" > "$(realpath PS-826-T.cnr)"
cnvkit.py segment --output PS-826-T.cns --processes 10 PS-826-T.cnr
cnvkit.py diagram --segment PS-826-T.cns --output Diagram.pdf PS-826-T.cnr && pdftoppm Diagram.pdf Diagram.png -png
cnvkit.py scatter --segment PS-826-T.cns --vcf "$(realpath ../week03/PS-826.PASS.vcf)" --output Scatter.pdf PS-826-T.cnr && pdftoppm Scatter.pdf Scatter -png
cnvkit.py heatmap --output Heatmap.pdf PS-826-T.cns && pdftoppm Heatmap.pdf Heatmap -png
