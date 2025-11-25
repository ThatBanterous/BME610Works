#!/bin/bash
[[ -e Plot ]] && rm -rf Plot
SigProfilerPlotting plotSBS --savefig_format pdf "$(realpath ./output/SBS/12.SBS96.exome)" "$(realpath Plot)" '12' '96'
SigProfilerPlotting plotSBS --savefig_format pdf "$(realpath ./output/SBS/12.SBS6.exome)" "$(realpath Plot)" '12' '6'
SigProfilerPlotting plotDBS --savefig_format pdf "$(realpath ./output/DBS/12.DBS78.exome)" "$(realpath Plot)" '12' '78'
