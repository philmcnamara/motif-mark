#! /home/phil/anaconda3/bin/python

# import cairo
import re
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
    '''All motifs for our search are stored in a list.'''
    motif_list = []
    with open(motifs_input, "r") as motifs_file:
        line = motifs_file.readline().strip()
        motif_list.append(line)
        for line in motifs_file:
            motif_list.append(line.strip())
    return motif_list


def parse_input_sequence(sequence_input):
    '''Reads input FASTA file into a list with the sequences in order.
    Each list index contains a dictionary with the key as the sequence
    ID from the header and the value as the base pair sequence.
    This preserves the order of input.'''
    sequence_list = []
    with open(sequence_input, "r") as sequence_file:
        line = sequence_file.readline().strip()
        while line:
            if(line[0] == ">"):
                sequence = []
                sequence.append(line[1:])
                sequence.append("")
                line = sequence_file.readline().strip()
                while line and line[0] != ">":
                    sequence[1] += line
                    line = sequence_file.readline().strip()
                sequence_list.append(sequence)
    return sequence_list


def find_motif_indices(motif_list, sequence_list):
    '''Finds provided motifs in the FASTA sequences with RegEx.
    Compatible with all IUPAC ambiguity codes.
    Returns a list of nested dictionaries, one dictionary per sequence ID.
    Each sequence ID dictionary has the sequence ID as a key
    and a dictionary of motifs as values.
    Motif values are a dictionary with the motif sequence as the key 
    and a list of indices as the values.
    These indices refer to left-most mapping positions of the 
    motif in the sequence.'''
    all_matches = []
    for sequence in sequence_list:
        sequence_ID = sequence[0]
        nucleotides = sequence[1]
        motif_matches = []
        any_matches = False
        for motif in motif_list:
            individual_match = False
            pattern = re.compile(motif)
            match_indices = []
            for match in pattern.finditer(nucleotides):
                any_matches = True
                individual_match = True
                match_indices.append(match.start())
            if individual_match:
                motif_matches.append([motif, match_indices])
        if any_matches:
            all_matches.append([sequence_ID, motif_matches])
    print(all_matches)


def draw_figures():
    pass


def main():
    args = get_arguments()
    motif_list = parse_input_motifs(args.motifs_input)
    sequence_list = parse_input_sequence(args.sequence_input)
    find_motif_indices(motif_list, sequence_list)

if __name__ == "__main__":
    main()
