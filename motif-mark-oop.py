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
        """Convert multi-line FASTA sequences into a dictionary with gene names as keys and 
        values as a list containing the sequence, start and end position, and strand."""
        if self.file_type == "fasta":
            with open(self.file, "r") as fin:
                sequence_dict = {}
                header = None
                seqs = []
                for line in fin:
                    line = line.strip()             
                    if line.startswith(">") == True and header == None:  # Header line
                        header = line.split()
                    elif line.startswith(">") != True:
                        seqs.append(line)
                    elif line.startswith(">") == True and header != None:
                        gene_name = header[0][1:]
                        location_full = header[1].split(":")
                        chrom = location_full[0]
                        location = location_full[1].split("-")
                        start_pos = location[0]
                        end_pos = location[1]
                        #join fasta sequence lines into one sequence
                        seqs = "".join(seqs)
                        #add fasta entry to dictionary
                        sequence_dict[gene_name] = [seqs, (chrom, start_pos, end_pos)]

                        #clear the seqs list after creating addition to dictionary
                        seqs = []
                        header = header = line.split()

                gene_name = header[0][1:]
                location_full = header[1].split(":")
                chrom = location_full[0]
                location = location_full[1].split("-")
                start_pos = location[0]
                end_pos = location[1]
                #join fasta sequence lines into one sequence
                seqs = "".join(seqs)
                #add fasta entry to dictionary
                sequence_dict[gene_name] = [seqs, (chrom, start_pos, end_pos)]
                #clear the seqs list after creating addition to dictionary
                seqs = []
                header = header = line.split()
            return(sequence_dict)

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
        self.seqtype = None
        self.introns = None
        self.exon = None
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
        self.exon = Exon(self.sequence[self.exon_start:self.exon_end])

    def intron_obj(self):
        '''generate intron objects from locations where exons are not present'''
        introns = []
        introns.append(self.sequence[0:self.exon_start])
        introns.append(self.sequence[self.exon_end:])
        self.introns =  introns

    def compare_motifs(self, motifs_list):
        '''Takes in a list of motifs and finds the locations of motifs in the ROI'''
        motif_positions = {}
        for motif_obj in motifs_list:
            sequence = self.sequence
            motif_positions[motif_obj.motif] = []
            for pos in range(len(sequence)):
                slide = sequence[pos:(pos+len(motif_obj.motif))]
                count = 0
                for i, locus in enumerate(slide):
                    if motif_obj.motif[i] == locus or motif_obj.motif[i].lower() == locus or motif_obj.motif[i].upper() == locus:
                        count += 1
                    elif motif_obj.motif[i] in motif_obj.nuc_notation_dict and locus in motif_obj.nuc_notation_dict[motif_obj.motif[i]]:
                        count += 1
                if count == len(motif_obj.motif):
                    motif_positions[motif_obj.motif].append((slide, (pos, pos+len(motif_obj.motif))))
        self.motif_positions = motif_positions

class Exon():
    def __init__(self, sequence):
        self.sequence = sequence       

class Motif():
    def __init__(self,motifseq):
        '''Creates a motif object'''
        self.motif = motifseq
        self.ambiguous_nts = None

    def nucleotide_notation(self):
        '''Generates a dictionary where keys = all possible ambiguous bases, and values = bases represented by that ambiguous base'''
        #Nucleic acid notation dictionary including DNA and RNA bases
        na_notation = { "W":("A","T","a","t"),"w":("A","T","a","t"),"S":("C","G","c","g"), "s":("C","G","c","g"),"Y":("C","T","c", "t"), "y":("c","t","C","T"),
        "M":("A","C","a","c"),"m":("A","C","a","c"),"K":("G","T","g","t"),"k":("G","T","g","t"), "R":("A","G","a","g"),"r":("A","G","a","g"), "B":("C","G","T","c","g","t"),"b":("C","G","T","c","g","t"),
        "D":("A","G","T","a","g","t"),"d":("A","G","T","a","g","t"), "H":("A","C","T","a","c","t"),"h":("A","C","T","a","c","t"), "V":("A","C","G","a","c","g"), "v":("A","C","G","a","c","g"), 
        "N":("A","C","G","T","a","c","g","t"),"n":("A","C","G","T","a","c","g","t"), 
        "U":("U","T","u","t"), "u":("U","T","u","t")}
        self.nuc_notation_dict = na_notation


class Drawing():
    def __init__(self, needed_for_drawing):
        self.to_draw = needed_for_drawing

    ##Methods##
    def draw_visualization(self):
        """Uses PyCairo to draw motifs on sequences and save as PNG."""
        number_drawings = len(self.to_draw)
        height_img = (number_drawings+1) * 200
        lengths = [len(seqs.sequence) for seqs in self.to_draw]
        longest_seq = max(lengths)
        width_img = longest_seq + 600
        WIDTH, HEIGHT = width_img, height_img
        
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
        ctx = cairo.Context(surface)
        ctx.set_source_rgb(1,1,1)
        ctx.paint()

        #start with lines for whole sequence (introns)
        for i, seq in enumerate(self.to_draw):
            sequence_len = len(seq.sequence)
            ctx.set_line_width(4)
            ctx.set_source_rgb(0,0,0)
            ctx.move_to(100,(200 * (i+1)))
            ctx.line_to(sequence_len + 100, (200 * (i+1)))
            ctx.stroke()
        
        #draw exons at the appropriate locations
        for i, seq in enumerate(self.to_draw):
            exon_start = seq.exon_start
            exon_end = seq.exon_end
            exonlength = exon_end - exon_start

            #exon box dimensions, colors, etc.
            ctx.set_source_rgb(0.18, 0.40, 0.64) #18, 40, 64
            ctx.rectangle(float(exon_start + 100), 175*(i+1)+(25*i), float(exonlength), 50)        
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
                    ctx.rectangle(float(motif_start + 100), 175 * (i+1)+(25*i), float(motif_end - motif_start), 50)        
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
            ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            ctx.move_to(100,150*(i+1)+(50*i))
            ctx.show_text(f"{gene_name} chr{chrom}:{start_pos} - {end_pos}")
            ctx.stroke()

        #make key for motifs present
        ctx.set_source_rgb(0.713, 0.886, 0.886) #182, 226, 226
        ctx.rectangle(float(longest_seq+200), 335, 300, 260)        
        ctx.fill()
        #title for legend
        ctx.set_source_rgb(0,0,0)
        ctx.set_font_size(28)
        ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        ctx.move_to(float(longest_seq+295),385)
        ctx.show_text("MOTIFS")
        #underline "MOTIFS" legend title
        ctx.set_line_width(2)
        ctx.set_source_rgb(0,0,0)
        ctx.move_to(float(longest_seq+295), 395)
        ctx.line_to(float(longest_seq+400), 395)
        ctx.stroke()

        #draw color coded motif boxes
        for i, seq in enumerate(self.to_draw):
            colors = [(0.5, 0, 0.3), (0.2, 0.7, 0.8), (1, 0.3, 0), (0.9, 0.7, 0.1), (0.2, 0.4, 0.8)]
            motifs = seq.motif_positions
            for j, motif in enumerate(motifs):

                #draw box for each motif within the key
                ctx.set_source_rgb(colors[j][0], colors[j][1], colors[j][2])
                ctx.rectangle(float(longest_seq+220), 365+((j+1)*45), float(30), 30)        
                ctx.fill()

                #draw text for label
                ctx.set_source_rgb(0,0,0)
                ctx.set_font_size(24)
                ctx.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
                ctx.move_to(float(longest_seq+280), 390+((j+1)*45))
                ctx.show_text(f"{motif.upper()}")      
            break

        #write to png output
        output_filename = fasta_file[:-6]
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