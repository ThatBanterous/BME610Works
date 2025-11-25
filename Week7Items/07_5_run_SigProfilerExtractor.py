#!/usr/bin/env python3
import argparse
from SigProfilerExtractor import sigpro


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input directory")
    parser.add_argument("output", help="Output directory")
    parser.add_argument("--cpu", help="Number of CPUs to use", type=int, default=5)

    return parser.parse_args()


def main():
    args = parse_arguments()

    sigpro.sigProfilerExtractor(input_type="vcf", input_data=args.input, output=args.output, reference_genome="GRCh38", exome=True, cpu=args.cpu, assignment_cpu=args.cpu)


if __name__ == "__main__":
    main()
