#!/usr/bin/env python3
import argparse
import os
import pandas
import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input MAF file")
    parser.add_argument("--gene", help="Select gene name", default="TP53")

    return parser.parse_args()


def main():
    args = parse_arguments()

    input_data = pandas.read_csv(args.input, sep="\t", comment="#")
    print(input_data)

    gene_data = input_data[(input_data["Hugo_Symbol"] == args.gene) & ~(input_data["HGVSp_Short"].isna())]
    print(gene_data)

    proteins = ""
    for _, row in tqdm.tqdm(gene_data.iterrows(), total=len(gene_data)):
        proteins += f"{row['HGVSp_Short'][2:]} "
    proteins = proteins.strip()

    with open(f"{args.gene}.sh", "w") as f:
        f.write(f"/BiO/Share/Tools/lollipops -legend -labels -o {args.gene}.png -dpi=600 -show-motifs {args.gene} {proteins}")

    os.system(f"bash {args.gene}.sh")


if __name__ == "__main__":
    main()
