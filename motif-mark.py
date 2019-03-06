#! /home/phil/anaconda3/bin/python

import cairo
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
    Returns their indices as 4 nested lists.
    Outermost list has all target sequences in order of input
    Each input sequence has its ID in index 0,
    and a list of motif matches in index 1.
    Each list of motif matches has the motif sequence in index 0 
    and a list of match positions in index 1
    '''
    all_matches = []
    for sequence in sequence_list:
        sequence_ID = sequence[0]
        nucleotides = sequence[1]
        motif_matches = []
        any_matches = False
        for motif in motif_list:
            individual_match = False
            pattern = re.compile(motif_to_regex(motif))
            match_indices = []
            for match in pattern.finditer(nucleotides):
                any_matches = True
                individual_match = True
                match_indices.append(match.start())
            if individual_match:
                motif_matches.append([motif, match_indices])
        if any_matches:
            all_matches.append([sequence_ID, motif_matches])
    return all_matches


def motif_to_regex(motif):
    '''Converts a motif with potential IUPAC ambiguity codes to a
    regex pattern which can match the desired sequence bases'''
    conversion_dict = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "U": "T",
        "M": "[AC]",
        "R": "[AG]",
        "W": "[AT]",
        "S": "[CG]",
        "Y": "[CT]",
        "K": "[GT]",
        "V": "[ACG]",
        "H": "[ACT]",
        "D": "[AGT]",
        "B": "[CGT]",
        "N": "[GATC]"
    }
    pattern = ""
    for base in motif:
        pattern += conversion_dict[base]
    return pattern


def draw_figures(matches):
    width = 1500
    # Length of matches = number of sequences, each sequence gets an independent figure
    height = 1000 * len(matches)
    surface = cairo.SVGSurface("figure.svg", width, height)
    ctx = cairo.Context(surface)

    # Draw sequence
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_line_width(10)
    ctx.move_to(100, 800)
    ctx.line_to(1400, 800)
    ctx.stroke()

    # Draw exon
    ctx.rectangle(600, 750, 300, 100) # x,y of top left corner, width, height
    ctx.fill()
    
    # Sequence ID Title
    ctx.move_to(100, 100)
    ctx.set_font_size(36)
    ctx.show_text("Sequence 1")

    # Legend
    ctx.move_to(1200, 100)
    ctx.show_text("Legend")

    # Divider to separate next figure
    ctx.move_to(0, 1000)
    ctx.set_source_rgb(0.5, 0.5, 0.5)
    ctx.line_to(1500, 1000)
    ctx.stroke()



def main():
    args = get_arguments()
    motif_list = ["GATAC", "MATAM"]
    sequence_list = parse_input_sequence(args.sequence_input)
    matches = find_motif_indices(motif_list, sequence_list)
    # print(matches[0]) # First sequence and matches
    # print(matches[0][0]) # First sequence ID
    # print(matches[0][1]) # First sequence list of matching motifs
    # print(matches[0][1][0]) # First sequence first motif, ID and indices
    # print(matches[0][1][0][0])  # First sequence first motif ID
    # print(matches[0][1][0][1])  # first sequence first motif list of indices
    draw_figures(matches)


if __name__ == "__main__":
    main()
