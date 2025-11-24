#!/usr/bin/env python3
"""
mutect.py: Mutect2 - Call somatic SNVs and indels via local assembly of haplotypes
"""
import argparse
import os
import sys
from utils import PipelineManagerBase
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))


class PipelineManager(PipelineManagerBase):
    def __init__(self, normal, tumor, output, config_file, dryrun):
        super().__init__(config_file, dryrun)
        self.normal = os.path.realpath(normal)
        self.tumor = os.path.realpath(tumor)
        self.output = os.path.realpath(output)
        self.normal_name = self.normal.split("/")[-1].split(".")[0]
        self.tumor_name = self.tumor.split("/")[-1].split(".")[0]
        self.name = self.output.split("/")[-1].split(".")[0]
        self.output_dir = os.path.dirname(self.output)

    def run_mutect(self, dependency_id=None):
        input_tumor = f"--input {self.tumor} --normal-sample {self.normal_name} "
        command = f"{self.config['TOOLS']['gatk']} Mutect2 --java-options \"{self.config['DEFAULT']['java_options']}\" --reference {self.config['REFERENCES']['fasta']} --input {self.normal} {input_tumor} --output {self.output_dir}/{self.name}.vcf --native-pair-hmm-threads {self.config['DEFAULT']['threads']} --max-mnp-distance 0
        self.create_sh("1.Mutect", command)
        return self.submit_job("1.Mutect", dependency_id=dependency_id)

    def run_filter(self, dependency_id=None):
        command = f"{self.config['TOOLS']['gatk']} FilterMutectCalls --java-options \"{self.config['DEFAULT']['java_options']}\" --reference {self.config['REFERENCES']['fasta']} --variant {self.output_dir}/{self.name}.vcf --output {self.output_dir}/{self.name}.filter.vcf"
        self.create_sh("2.Filter", command)
        return self.submit_job("2.Filter", dependency_id=dependency_id, cpus=1)

    def run_pass(self, dependency_id=None):
        command = f"{self.config['TOOLS']['awk']} -F '\t' '{{if($0 ~ /\\#/) print; else if($7 == \"PASS\") print}}' {self.output_dir}/{self.name}.filter.vcf > {self.output_dir}/{self.name}.PASS.vcf"
        self.create_sh("3.PASS", command)
        return self.submit_job("3.PASS", dependency_id=dependency_id, cpus=1)

    def run_index(self, dependency_id=None):
        command = f"{self.config['TOOLS']['gatk']} IndexFeatureFile --java-options \"{self.config['DEFAULT']['java_options']}\" --input {self.output_dir}/{self.name}.PASS.vcf --output {self.output_dir}/{self.name}.PASS.vcf.idx"
        self.create_sh("4.Index", command)
        return self.submit_job("4.Index", dependency_id=dependency_id, cpus=1)

    def run_maf(self, dependency_id=None):
        command = f"{self.config['TOOLS']['vcf2maf']} --vep-path {self.config['TOOLS']['vep']} --vep-data {self.config['TOOLS']['vep']} --vep-forks {self.config['DEFAULT']['threads']} --ncbi-build 'GRCh38' --input-vcf {self.output_dir}/{self.name}.PASS.vcf --output {self.output_dir}/{self.name}.PASS.maf --tumor-id {self.tumor_name} --normal-id {self.normal_name} --ref-fasta {self.config['REFERENCES']['fasta']} --vep-overwrite"
        self.create_sh("5.MAF", command)
        return self.submit_job("5.MAF", dependency_id=dependency_id)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("normal", help="Normal BAM file")
    parser.add_argument("tumor", help="Tumor BAM file")
    parser.add_argument("output", help="Output MAF file")
    parser.add_argument("-c", "--config", help="config INI file", default="config.ini")
    parser.add_argument("-n", "--dryrun", help="Don't actually run any recipe; just make .SH only", default=False, action="store_true")

    return parser.parse_args()


def main():
    args = parse_arguments()

    pipeline = PipelineManager(normal=args.normal, tumor=args.tumor, output=args.output, config_file=args.config, dryrun=args.dryrun)

    pipeline.create_dir()

    mutect_job_id = pipeline.run_mutect()
    filter_job_id = pipeline.run_filter(dependency_id=mutect_job_id)
    pass_job_id = pipeline.run_pass(dependency_id=filter_job_id)
    index_job_id = pipeline.run_index(dependency_id=pass_job_id)
    pipeline.run_maf(dependency_id=index_job_id)


if __name__ == "__main__":
    main()
