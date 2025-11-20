#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610-04
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
python3 3_draw_rainfall.py "$(realpath ../SAMPLE.PASS.vcf)" "$(realpath .)"
