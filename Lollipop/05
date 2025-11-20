#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610-05
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
python3 3_lollipop.py "$(realpath ../SAMPLE.PASS.maf)" --gene 'TP53'
