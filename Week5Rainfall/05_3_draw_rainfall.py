#!/usr/bin/env python3
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input merged VCF file")
    parser.add_argument("output", help="Output directory")

    return parser.parse_args()


def reverse_strand(code):
    match_nucleotide = {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">"}
    return "".join(list(map(lambda x: match_nucleotide[x], code)))


def main():
    args = parse_arguments()

    parameters = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 75, "xtick.labelsize": 50, "ytick.labelsize": 50, "legend.fontsize": 30, "legend.title_fontsize": 30, "figure.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42}
    matplotlib.use("Agg")
    matplotlib.rcParams.update(parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=parameters)

    names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Normal", "Tumor"]
    input_data = pandas.read_csv(args.input, sep="\t", header=None, comment="#", names=names)
    print(input_data)

    mutation_set = {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"}
    mutation_list = []
    for _, row in tqdm.tqdm(input_data.iterrows(), total=len(input_data)):
        ref = row["REF"]
        alt = row["ALT"]
        mutation = f"{ref}>{alt}"

        if (ref == "G") or (ref == "A"):
            mutation = reverse_strand(mutation)

        if mutation not in mutation_set:
            mutation_list.append("Indel")
        else:
            mutation_list.append(mutation)
    input_data["Mut"] = mutation_list
    print(input_data)

    previous_position = 999999999
    distance_list = []
    for _, row in tqdm.tqdm(input_data.iterrows(), total=len(input_data)):
        distance = row["POS"] - previous_position

        if distance < 0:
            distance_list.append(0)
        else:
            distance_list.append(distance)

        previous_position = row["POS"]
    input_data["Dis"] = distance_list
    print(input_data)

    chromosome_list = sorted(set(input_data["CHROM"]))
    print(len(chromosome_list), chromosome_list)

    coloring = {"C>A": "tab:blue", "C>G": "tab:orange", "C>T": "tab:green", "T>A": "tab:purple", "T>C": "tab:brown", "T>G": "tab:olive", "Indel": "tab:gray"}

    for chromosome in tqdm.tqdm(chromosome_list):
        chromosome_data = input_data.loc[(input_data["CHROM"] == chromosome) & (input_data["Dis"] != 0)]

        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        seaborn.scatterplot(data=chromosome_data, x="POS", y="Dis", hue="Mut", style="Mut", palette=coloring, legend="full", s=1000, ax=ax)

        matplotlib.pyplot.title(chromosome)
        matplotlib.pyplot.xlabel("Mutation position (bp)")
        matplotlib.pyplot.ylabel("Distance (bp)")
        matplotlib.pyplot.yscale("log", base=10)

        matplotlib.pyplot.tight_layout()
        fig.savefig(f"{chromosome}.png")
        matplotlib.pyplot.close(fig)


if __name__ == "__main__":
    main()
