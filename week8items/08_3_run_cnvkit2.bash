#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=week08
#SBATCH --mem=20G
#SBATCH --export=ALL

# --- Step 1: Run the Batch Command ---
# We use your 12N (Normal) and 12T (Tumor) BAM files.
# Note: The output files will automatically be named after the Tumor file (12T.Sort.MarkDuplicates.BQSR...)
cnvkit.py batch \
    --processes 10 \
    --normal "$(realpath 12N.Sort.MarkDuplicates.BQSR.bam)" \
    --fasta "/BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta" \
    --targets "/BiO/Share/Database/UCSC/hg38/UCSC_hg38.bed" \
    --annotate "/BiO/Share/Database/UCSC/hg38/refFlat.txt"  \
    --scatter \
    --diagram \
    --output-dir "$(realpath .)" \
    "$(realpath 12T.Sort.MarkDuplicates.BQSR.bam)"

# --- Step 2: Clean the Data ---
# The batch command outputs a file named "12T.Sort.MarkDuplicates.BQSR.cnr".
# We filter out non-standard chromosomes (those with '_' in the name) and save it as a simpler name: "12T.cnr"
grep --invert-match '_' "$(realpath ./12T.Sort.MarkDuplicates.BQSR.cnr)" > "$(realpath 12T.cnr)"

# --- Step 3: Segmentation ---
# We use the clean "12T.cnr" file
cnvkit.py segment --output 12T.cns --processes 10 12T.cnr

# --- Step 4: Plotting ---
# 1. Diagram (Whole Genome View)
cnvkit.py diagram --segment 12T.cns --output Diagram.pdf 12T.cnr && pdftoppm Diagram.pdf Diagram -png

# 2. Scatter Plot (Detailed View with VCF overlay)
# Note: We point to your 12.PASS.vcf from the earlier step
cnvkit.py scatter --segment 12T.cns --vcf "$(realpath 12.PASS.vcf)" --output Scatter.pdf 12T.cnr && pdftoppm Scatter.pdf Scatter -png

# 3. Heatmap
cnvkit.py heatmap --output Heatmap.pdf 12T.cns && pdftoppm Heatmap.pdf Heatmap -png
