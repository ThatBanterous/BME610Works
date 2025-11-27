#!/bin/bash
pip3 install --requirement "$(realpath 08_2_requirements.txt)"

# FIXED: Added 'repos' to specify the download server (mirror) and 'ask=FALSE' to prevent hanging.
Rscript --vanilla -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")' -e 'BiocManager::install("DNAcopy", update=FALSE, ask=FALSE)'
