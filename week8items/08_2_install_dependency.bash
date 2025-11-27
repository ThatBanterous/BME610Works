#!/bin/bash
pip3 install --requirement "$(realpath 08_2_requirements.txt)"
Rscript --vanilla -e 'install.packages("BiocManager")' -e 'BiocManager::install("DNAcopy")'
