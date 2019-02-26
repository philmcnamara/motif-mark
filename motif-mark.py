#! /home/phil/anaconda3/bin/python

# import cairo
# import re
import argparse
# import math


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Draws a sequence region with motifs annotated"
    )
    # Motif Input
    parser.add_argument(
        "-m", "--motifs_input",
        help="input file with one motif per line, max 10",
        required=True, type=str
    )
    # FASTA file with sequences
    parser.add_argument(
        "-s", "--sequence_input",
        help="FASTA input file with sequences, max 10",
        required=True, type=str
    )
    return parser.parse_args()


def parse_input_motifs(motifs_input):
    motif_list = []
    with open(motifs_input, "r") as motifs_file:
        line = motifs_file.readline().strip()
        motif_list.append(line)
        for line in motifs_file:
            motif_list.append(line.strip())
    return motif_list


def parse_input_sequence(sequence_input):
    sequence_list = []
    with open(sequence_input, "r") as sequence_file:
        line = sequence_file.readline().strip()
        while line:
            if(line[0] == ">"):
                sequence = {}
                sequence_name = line[1:]
                sequence[sequence_name] = ""
                sequence_list.append(sequence)
            line = sequence_file.readline().strip()
    return sequence_list


def find_motif_indices():
    pass


def draw_figures():
    pass


def main():
    args = get_arguments()
    print(parse_input_motifs(args.motifs_input))
    print(parse_input_sequence(args.sequence_input))


if __name__ == "__main__":
    main()
