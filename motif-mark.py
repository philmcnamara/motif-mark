#! /home/phil/anaconda3/bin/python

import cairo
import re
import argparse


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
        help="FASTA input file with sequences",
        required=True, type=str
    )
    # Output figure name
    parser.add_argument(
        "-o", "--output",
        help="Output file name (omit extension)",
        required=False, default="Figure"
    )
    # Scale bar
    parser.add_argument(
        "-l", "--scale_bar_length",
        help="length of the scale bar in bp (base pairs)",
        required=False, default=100, type=int
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
            pattern = re.compile(motif_to_regex(motif.upper()))
            match_indices = []
            for match in pattern.finditer(nucleotides.upper()):
                any_matches = True
                individual_match = True
                match_indices.append(match.start())
            if individual_match:
                motif_matches.append([motif, match_indices])
        if any_matches:
            all_matches.append([sequence_ID, nucleotides, motif_matches])
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


def draw_figures(matches, output, scale_bar_length):
    width = 1500
    # 1000 pixels per subfigure, 1 per sequence
    height = 1000 * len(matches)
    surface = cairo.SVGSurface(output + ".svg", width, height)
    ctx = cairo.Context(surface)

    # RGB colors for marking the motifs
    colors = [
        [0.2, 0.7, 0.2],   # green
        [0.7, 0.45, 0.2],  # orange
        [0.57, 0.2, 0.7],  # purple
        [0.7, 0.7, 0.2],   # yellow
        [0.2, 0.2, 0.7],   # dark blue
        [0.7, 0.2, 0.2],   # red
        [0.7, 0.7, 0.7],   # lavender
        [0.2, 0.7, 0.7],   # sky blue
        [0.2, 0.7, 0.5],   # teal
        [0.5, 0.5, 0.5]    # gray
    ]
    
    for s in range(len(matches)):
        sequence = matches[s][1]

        # Sequence is always drawn with a length of 1300 pixels
        # scale_factor is pixels/base
        scale_factor = 1300/len(sequence)

        # Draw sequence
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(10)
        ctx.move_to(100, 800 + 1000 * s)
        ctx.line_to(1400, 800 + 1000 * s)
        ctx.stroke()

        # Find exon start and end position
        start = 0
        end = 0
        for i in range(len(sequence)):
            if sequence[i].isupper():
                start = i
                break
        # Iterate backwards to find the end
        for i in range(len(sequence) - 1, 0, -1):
            if sequence[i].isupper():
                end = i
                break

        # Draw exon
        # x,y of top left corner, width, height
        # offset by 100 for margins
        ctx.rectangle(start*scale_factor + 100, 775 + 1000 * s,
                      (end-start)*scale_factor, 50)
        ctx.fill()

        # Draw motifs
        motifs = matches[s][2]
        for i in range(len(motifs)):
            # Each new motif gets the next set of RBG values from colors
            ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            for match_position in range(len(motifs[i][1])):
                start_position = motifs[i][1][match_position]
                # offset by 100 for margins
                ctx.rectangle(start_position*scale_factor + 100,
                              775 + 1000 * s, 10, 50)
                ctx.fill()

        # Sequence ID Title
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(100, 100 + 1000 * s)
        ctx.set_font_size(36)
        ctx.show_text(matches[s][0])

        # Legend
        ctx.move_to(1200, 100 + 1000 * s)
        ctx.show_text("Motif Legend")
        ctx.set_font_size(24)
        for i in range(len(motifs)):
            ctx.move_to(1250, 160 + i * 50 + 1000 * s)
            ctx.show_text(motifs[i][0])
            ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            ctx.rectangle(1210, 135 + i * 50 + 1000 * s, 25, 25)
            ctx.fill()
            ctx.set_source_rgb(0, 0, 0)

        # Scale Bar
        ctx.move_to(100, 275 + 1000 * s)
        ctx.set_font_size(36)
        ctx.show_text("Scale")
        ctx.set_font_size(24)
        ctx.move_to(100, 325 + 1000 * s)
        ctx.show_text(str(scale_bar_length) + "bp")
        ctx.move_to(100, 350 + 1000 * s)
        ctx.line_to(100 + scale_factor * scale_bar_length, 350 + 1000 * s)
        ctx.stroke()

        # Divider to separate next figure
        ctx.move_to(0, 1000 + 1000 * s)
        ctx.set_source_rgb(0.5, 0.5, 0.5)
        ctx.line_to(1500, 1000 + 1000 * s)
        ctx.stroke()


def main():
    args = get_arguments()
    sequence_list = parse_input_sequence(args.sequence_input)
    motif_list = parse_input_motifs(args.motifs_input)
    matches = find_motif_indices(motif_list, sequence_list)
    draw_figures(matches, args.output, args.scale_bar_length)


if __name__ == "__main__":
    main()
