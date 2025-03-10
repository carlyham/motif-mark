#!/usr/bin/env python

import argparse
import cairo

def get_args():
    parser = argparse.ArgumentParser(description="Motif Marking Visualization")
    parser.add_argument('-f', type=str, help='Input FASTA file', required = True)
    parser.add_argument('-m', type=str, help='Input motifs file', required = True)
    return parser.parse_args()
args = get_args()

fasta_file = args.f
motif_file = args.m

class ParseInputFile:
    def __init__(self, f_file):
        """Creates objects of the input file class, motif or fasta"""
        self.file = f_file
        self.file_type = None       

    ## Methods ##
    def what_file_type(self):
        '''Determines what type of file is being used as input, fasta or motif'''
        if self.file.endswith(".fasta") == True:
            self.file_type = "fasta"
        elif self.file.endswith(".txt") == True:
            self.file_type = "motif"

    def parse_fasta_file(self):
        """Convert multi-line FASTA sequences into a dictionary with gene names as keys 
    and values as a list containing the sequence, start and end position, and strand."""
        if self.file_type == "fasta":
            with open(self.file, "r") as fin:
                sequence_dict = {}
                header = None
                seqs = []

                for line in fin:
                    line = line.strip() 
                    # process new header line
                    if line.startswith(">"):
                        if header is not None:
                            gene_name = header[0][1:]
                            location_full = header[1].split(":")
                            chrom = location_full[0]
                            location = location_full[1].split("-")
                            start_pos = location[0]
                            end_pos = location[1]
                            # add seq lines into a single sequence
                            seqs = "".join(seqs)
                            sequence_dict[gene_name] = [seqs, (chrom, start_pos, end_pos)]
                            # empty seqs to continue looping
                            seqs = []
                        # add new header to header variable
                        header = line.split()
                    # add seq lines to seqs list
                    else:
                        seqs.append(line)
                # process the last entry after the loop
                if header is not None:
                    gene_name = header[0][1:]
                    location_full = header[1].split(":")
                    chrom = location_full[0]
                    location = location_full[1].split("-")
                    start_pos = location[0]
                    end_pos = location[1]
                    seqs = "".join(seqs)
                    sequence_dict[gene_name] = [seqs, (chrom, start_pos, end_pos)]
        return sequence_dict


    def parse_motif_file(self):
        '''Take text file containing motifs into a list'''
        motifs = []
        if self.file_type == "motif":
            with open(self.file, "r") as mf:
                for line in mf:
                    line = line.strip()
                    motifs.append(line)
        return motifs



class RegionOfInterest():
    def __init__(self, sequence, gene_name, chrom, start_pos, end_pos):
        '''Creates an object of class region of interest to find introns and exons, align motifs.'''
        self.sequence = sequence
        self.gene_name = gene_name
        self.chrom = chrom
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.type_of_seq = None
        self.introns = None
        self.exons = None
        self.motif_positions = None

    ##Methods##
    def find_exons(self):
        '''Get exon start position for pycairo drawings'''
        exon_start = None
        exon_end = None
        for pos, base in enumerate(self.sequence):
            if base.isupper() == True and exon_start == None:
                exon_start = pos + 1
            elif exon_start != None and exon_end == None and base.islower() == True:
                exon_end = pos + 1
        self.exon_start = exon_start
        self.exon_end = exon_end
    def exon_obj(self):
        '''generate a new object of the exon class to save as an attribute attached to region of interest'''
        self.exons = Exon(self.sequence[self.exon_start:self.exon_end])

    def intron_obj(self):
        '''generate intron objects from locations where exons are not present'''
        introns = []
        introns.append(self.sequence[0:self.exon_start])
        introns.append(self.sequence[self.exon_end:])
        self.introns =  introns

    def compare_motifs(self, motifs_list):
        '''Takes a list of motifs and finds the locations of motifs in the sequence'''
        motif_positions = {}
        # Loop through each motif object
        for motif_obj in motifs_list:
            motif = motif_obj.motif
            motif_positions[motif] = []
            
            # Slide through the sequence
            for pos in range(len(self.sequence) - len(motif) + 1):
                slide = self.sequence[pos:(pos+len(motif))]
                is_match = True

                # Compare each base
                for i in range(len(motif)):
                    motif_base = motif[i]
                    seq_base = slide[i]
                    # Handle exact or case-insensitive matches
                    if motif_base.upper() == seq_base.upper():
                        continue
                    # Handle ambiguous nucleotide matches
                    elif motif_base in motif_obj.nuc_notation_dict:
                        if seq_base not in motif_obj.nuc_notation_dict[motif_base]:
                            is_match = False
                            break
                    else:
                        is_match = False
                        break
                # If a full match is found, add it to the list
                if is_match:
                    motif_positions[motif].append((slide, (pos, pos+len(motif))))
        self.motif_positions = motif_positions


class Exon():
    def __init__(self, sequence):
        self.sequence = sequence       


class Motif():
    def __init__(self,motifseq):
        '''Creates a motif object'''
        self.motif = motifseq

    def nucleotide_notation(self):
        '''Generates a dictionary where keys = all possible ambiguous bases, and values = bases represented by that ambiguous base'''
        #Nucleic acid notation dictionary including DNA and RNA bases
        na_notation = { "W":("A","T","a","t"),"w":("A","T","a","t"), "S":("C","G","c","g"), "s":("C","G","c","g"),"M":("A","C","a","c"),"m":("A","C","a","c"),
        "K":("G","T","g","t"),"k":("G","T","g","t"), "N":("A","C","G","T","a","c","g","t"),"n":("A","C","G","T","a","c","g","t"),"R":("A","G","a","g"),"r":("A","G","a","g"),
        "Y":("C","T","c", "t"), "y":("c","t","C","T"), "B":("C","G","T","c","g","t"), "b":("C","G","T","c","g","t"),
        "D":("A","G","T","a","g","t"),"d":("A","G","T","a","g","t"), "H":("A","C","T","a","c","t"),"h":("A","C","T","a","c","t"), 
        "V":("A","C","G","a","c","g"), "v":("A","C","G","a","c","g"), "U":("U","T","u","t"), "u":("U","T","u","t")}
        self.nuc_notation_dict = na_notation


class Drawing():
    def __init__(self, needed_for_drawing):
        self.to_draw = needed_for_drawing

    ##Methods##
    def draw_visualization(self):
        """Uses PyCairo to draw motifs on sequences and save as PNG.""" 
        lengths = [len(seqs.sequence) for seqs in self.to_draw]
        longest_seq = max(lengths)
        width_img = longest_seq + 600
        number_drawings = len(self.to_draw)
        height_img = (number_drawings+1) * 180
        WIDTH, HEIGHT = width_img, height_img
        
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
        ctx = cairo.Context(surface)
        ctx.set_source_rgb(1,1,1)
        ctx.paint()

        #start with lines for whole sequence (introns)
        for i, seq in enumerate(self.to_draw):
            sequence_len = len(seq.sequence)
            ctx.set_line_width(4)
            ctx.set_line_cap(cairo.LINE_CAP_ROUND)
            ctx.set_source_rgb(0,0,0)
            ctx.move_to(100,(200 * (i+1)))
            ctx.line_to(sequence_len + 100, (200 * (i+1)))
            ctx.stroke()
        
        #draw exons at the appropriate locations
        for i, seq in enumerate(self.to_draw):
            exon_start = seq.exon_start
            exon_end = seq.exon_end
            exonlength = exon_end - exon_start

            #exons box dimensions, colors, etc.
            ctx.set_source_rgb(0.18, 0.40, 0.64) #18, 40, 64
            ctx.rectangle((exon_start + 100), 174*(i+1)+(28*i), (exonlength), 40)        
            ctx.fill()
        
        #draw motifs, color code by motif
        for i, seq in enumerate(self.to_draw):
            #color code key
            colors = [(0.5, 0, 0.3), (0.2, 0.7, 0.8), (1, 0.3, 0), (0.9, 0.7, 0.1), (0.2, 0.4, 0.8)]
            motifs = seq.motif_positions
            #iterate through motifs
            for j, motif in enumerate(motifs):
                motifcolor = colors[j]
                ctx.set_source_rgb(motifcolor[0], motifcolor[1], motifcolor[2])
                for location in motifs[motif]:
                    motif_start = location[1][0]
                    motif_end = location[1][1]
                    ctx.rectangle((motif_start + 100), 174 * (i+1)+(28*i), (motif_end - motif_start), 40)        
                    ctx.fill()

        #label each sequence with gene name, start and end position, and chromosome
        for i, seq in enumerate(self.to_draw):
            gene_name = seq.gene_name
            chrom = seq.chrom
            chrom = chrom.split("r")
            chrom = chrom[1]
            start_pos = seq.start_pos
            end_pos = seq.end_pos

            ctx.set_source_rgb(0, 0, 0)
            ctx.set_font_size(25)
            ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL)
            ctx.move_to(100,150*(i+1)+(50*i))
            ctx.show_text(f"{gene_name} chr{chrom}:{start_pos} - {end_pos}")
            ctx.stroke()

        output_filename = fasta_file[:-6]
        #title for motif mark image
        ctx.set_source_rgb(0,0,0)
        ctx.set_font_size(35)
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL)
        ctx.move_to(75, 75)
        ctx.show_text(f"Motifs Present in {output_filename} Genes of Interest")

        #title for legend
        ctx.set_source_rgb(0,0,0)
        ctx.set_font_size(28)
        ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL)
        ctx.move_to((longest_seq+295),385)
        ctx.show_text("MOTIFS")
        #underline "MOTIFS" legend title
        ctx.set_line_width(2)
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.set_source_rgb(0,0,0)
        ctx.move_to((longest_seq+295), 395)
        ctx.line_to((longest_seq+400), 395)
        ctx.stroke()

        #draw color coded motif boxes
        for i, seq in enumerate(self.to_draw):
            colors = [(0.5, 0, 0.3), (0.2, 0.7, 0.8), (1, 0.3, 0), (0.9, 0.7, 0.1), (0.2, 0.4, 0.8)]
            motifs = seq.motif_positions
            for j, motif in enumerate(motifs):

                #draw box for each motif within the key
                ctx.set_source_rgb(colors[j][0], colors[j][1], colors[j][2])
                ctx.rectangle((longest_seq+220), 365+((j+1)*45), 30, 30)        
                ctx.fill()

                #draw text for label
                ctx.set_source_rgb(0,0,0)
                ctx.set_font_size(24)
                ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
                ctx.move_to((longest_seq+280), 390+((j+1)*45))
                ctx.show_text(f"{motif.upper()}")      
            break

        #write to png output
        surface.write_to_png(f"{output_filename}.png")


##Using the actual inputs to create objects
input_fasta = ParseInputFile(fasta_file)
input_fasta.what_file_type()
fasta_dict = input_fasta.parse_fasta_file()


input_motif = ParseInputFile(motif_file)
input_motif.what_file_type()
motif_list = input_motif.parse_motif_file()


seqs_obj = [RegionOfInterest(fasta_dict[seqs][0], seqs, fasta_dict[seqs][1][0], fasta_dict[seqs][1][1], fasta_dict[seqs][1][2]) for seqs in fasta_dict]


motifs_lst = [Motif(x) for x in motif_list]

for motif_obj in motifs_lst:
    motif_obj.nucleotide_notation()


for sequence_obj in seqs_obj:
    sequence_obj.find_exons()
    sequence_obj.exon_obj()
    sequence_obj.intron_obj()
    sequence_obj.compare_motifs(motifs_lst)

final_drawing = Drawing(seqs_obj)
final_drawing.draw_visualization()