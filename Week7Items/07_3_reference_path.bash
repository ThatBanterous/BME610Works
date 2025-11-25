#!/bin/bash
[[ ! -e "$(realpath ./lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/matrix)" ]] && ln -sv /BiO/Live/jwlee230/BME610/week07/lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/matrix "$(realpath ./lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/matrix)"
[[ ! -e "$(realpath ./lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/vcf_files)" ]] && ln -sv /BiO/Live/jwlee230/BME610/week07/lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/vcf_files "$(realpath ./lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/vcf_files)"

[[ -e "$(realpath ./input)" ]] && rm -rfv "$(realpath ./input)"
mkdir input && ln -sv "$(realpath 12.PASS.vcf)" "$(realpath input)"
