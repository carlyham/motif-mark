#!/usr/bin/env python

import argparse
import cairo

def get_args():
    parser = argparse.ArgumentParser(description="Motif Marking Visualization")
    parser.add_argument('-f', type=str, help='Input FASTA file', default='/Users/carlyhamilton/bioinfo/Bi625/motif_mark/Figure_1.fasta')
    parser.add_argument('-m', type=str, help='Input motifs file', default='/Users/carlyhamilton/bioinfo/Bi625/motif_mark/Fig_1_motifs.txt')
    return parser.parse_args()
args = get_args()

fasta_file = args.f
motif_file = args.m

class InputFastaFile:
    def __init__(self, fasta_file):
        """Initialize with a FASTA file and store sequences as a dictionary."""
        self.file = fasta_file
        self.sequence = {}  # Store sequences
        

    ## Methods ##
    def oneline_fasta_list(self):
        """Convert multi-line FASTA sequences into a dictionary with headers as keys."""
        sequence = {}  # Temporary dictionary
        header = None
        seq = ""     
        with open(self.file, "r") as fin:
            for line in fin:
                line = line.strip()  # Remove whitespace             
                if line.startswith(">"):  # Header line
                    if header:  # Store the previous sequence before moving to a new one
                        sequence[header] = seq
                    header = line  # Update header
                    seq = ""  # Reset sequence storage
                else:
                    seq += line  # Append sequence data 
        # Save the last sequence after exiting the loop
        if header:
            sequence[header] = seq 
        self.sequence = sequence  # Store the dictionary in the object
        return sequence  # Return the dictionary

    

class InputMotifsFile:
    def __init__(self, motif_file):
        '''Takes in a txt file of motifs and returns a list of all possible motifs.'''
        ##Data##
        self.file = motif_file
        self.expanded_sequences = []
        

    ##Methods##
    def parse_motifs_file(self, file):
        motifs = []
        with open(file, "r") as fh:
            for line in fh:
                line = line.strip()
                line = line.upper()
                motifs.append(line)
        return motifs

    def expand_motif(seq, na_notation, index=0, current=""):
        """Recursively expands a sequence with ambiguous bases."""
        if index == len(seq):  # Do this when full sequence generated, returns current seq as list 
            return [current]

        char = seq[index]
        expanded_sequences = []

        if char in na_notation:  # If ambiguous, replace with possible values
            for replacement in na_notation[char]:
                expanded_sequences.extend(expand_motif(seq, na_notation, index + 1, current + replacement))
        else:  # If not ambiguous, keep the character
            expanded_sequences.extend(expand_motif(seq, na_notation, index + 1, current + char))

        self.expanded_sequences = expanded_sequences
        return expanded_sequences

#Nucleic acid notation dictionary including DNA and RNA bases
na_notation = {"W":["W","A","T"], "S":["S","C","G"], "M":["M","A", "C"], "K":["K","G","T"],
            "N": ["N","A","G","T","C"], "R": ["R","A","G"] , "Y":["Y","C","T"], "U":["T"],
            "B":["B","C", "G", "T"], "D": ["D","A", "T", "G"], "H": ["H","A","T","C"], "V": ["V","A","G","C"]}


            

class RegionOfInterest(InputFastaFile, InputMotifsFile):
    def __init__(self, expanded_sequences): 
        InputFastaFile.__init__(sequence)
        InputMotifsFile.__init__(expanded_sequences)
        self.length = 0
        self.intron_or_exon = []
        self.motif_locations = {}

    def find_introns_exons(self, sequence):
        intron_or_exon = []
        self.length = len(sequence)
        for pos in range(1, len(sequence)):
            if sequence[pos-1].isupper() and sequence[pos].islower():
                intron_or_exon.append((pos,"exon"))
            elif sequence[pos-1].islower() and sequence[pos].isupper():
                intron_or_exon.append((pos,"intron"))
        #list of locations where sequence switches FROM being an intron to exon and vice versa.
        self.intron_or_exon = intron_or_exon
        return intron_or_exon

    def compare_motifs(self, sequence, expanded_sequences):
        motif_positions = {}
        for motif in expanded_sequences:
            motif_size = len(motif)
            for i, pos in enumerate(sequence):
                if motif in sequence and (motif, motif_size) not in motif_positions.keys:
                    motif_positions[(motif, motif_size)] = [i]
                elif motif in sequence and (motif, motif_size) in motif_positions.keys:
                    motif_positions[(motif, motif_size)].append(i)
        self.motif_locations = motif_positions
        return motif_positions



class Drawing(RegionOfInterest):
    def __init__(self, fasta_file, motif_file):
        super().__init__(fasta_file, motif_file)

    def draw_visualization(self, output_file):
        """Uses PyCairo to draw motifs on sequences and save as PNG."""
        width, height = 1000, (len(self.sequence) * 100) + 100
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)
        
        ctx.set_font_size(16)
        y_offset = 50

        colors = [(0.5, 0, 0.3), (0.2, 0.7, 0.8), (1, 0.3, 0), (0.9, 0.7, 0.1), (0.2, 0.4, 0.8)]  # Different motif colors

        for i, (header, seq) in enumerate(self.sequence.items()):
            ctx.set_source_rgb(0, 0, 0)  # Black for sequence line
            ctx.move_to(50, y_offset + (i * 100))
            ctx.line_to(950, y_offset + (i * 100))
            ctx.stroke()

            motif_positions = self.compare_motifs(seq, self.motifs)

            for idx, (motif, positions) in enumerate(motif_positions.items()):
                ctx.set_source_rgb(*colors[idx % len(colors)])
                for start, end in positions:
                    box_start = 50 + (900 * start / len(seq))
                    box_end = 50 + (900 * end / len(seq))
                    ctx.rectangle(box_start, y_offset + (i * 100) - 5, box_end - box_start, 10)
                    ctx.fill()

            ctx.set_source_rgb(0, 0, 0)
            ctx.move_to(10, y_offset + (i * 100))
            ctx.show_text(header)

        surface.write_to_png(output_file)
        print(f"Visualization saved as {output_file}")


#Read and process the FASTA file
fasta_obj = InputFastaFile(args.f)
sequences = fasta_obj.oneline_fasta_list()  # Dictionary of {header: sequence}

#Read and expand motifs
motif_obj = InputMotifsFile(args.m)
all_motifs = motif_obj.parse_motifs_file()

# Expand ambiguous motifs
expanded_motifs = []
for motif in all_motifs:
    expanded_motifs.extend(motif_obj.expand_motifs(motif))

#find introns/exons and motif positions
seqs_objects = {}  # Store results per sequence

for header, sequence in sequences.items():
    roi = RegionOfInterest()
    roi.find_introns_exons(sequence)
    motif_positions = roi.compare_motifs(sequence, expanded_motifs)
    seqs_objects[header] = roi  # Store results

# Step 4: Draw the visualization
drawing = Drawing(expanded_motifs)  
drawing.sequences = sequences  # Attach sequences for drawing
drawing.motifs = expanded_motifs  # Attach motifs for visualization
output_file = fasta_file.rsplit(".", 1)[0] + ".png"  # Match input prefix
drawing.draw_visualization(output_file)