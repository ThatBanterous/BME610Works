#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

python3 -B mutect.py $(realpath 12N.Sort.MarkDuplicates.BQSR.bam) $(realpath 12T.Sort.MarkDuplicates.BQSR.bam) 12.maf
