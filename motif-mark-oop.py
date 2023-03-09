#!/usr/bin/env python3.10
import argparse
import re
import cairo

def get_args():
    parser = argparse.ArgumentParser(description = "Program to visualize motifs and exons within fasta sequence.")
    parser.add_argument("-f","--fasta", help ="Fasta file to search.", type = str, required = True)
    parser.add_argument("-m","--motif",help = "File containing target motifs.",type = str, required = True)
    return parser.parse_args()
args = get_args()

fasta_file = args.fasta
motif_file = args.motif

import re

class FastaParse:
    def __init__(self, filename):
        self.filename = filename
        self.records = []

    def parse_fasta(self):
        'Parses a fasta file to separate sequence and header lines. Storing them in a list of tuples. Concatenates seq lines into one'
        with open(self.filename) as file:
            record = {'header': '', 'sequence': '', 'exons': []}
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if record['header']:
                        self.records.append(record)
                        record = {'header': '', 'sequence': '', 'exons': []}
                    record['header'] = line[1:]
                else:
                    record['sequence'] += line
                    if 'exon' in record['header'].lower():
                        # Extract start and end positions of exon from the header
                        exon_start, exon_end = re.findall(r'\d+', record['header'])
                        # Add exon positions to the list of exons for this record
                        record['exons'].append((int(exon_start), int(exon_end)))
            if record['header']:
                self.records.append(record)
        return self.records

    def header_info(self):
        'Used to extract individual components of fasta headers. Name, start, stop. Can be modified as needed.'
        for record in self.records:
            header = record['header']
            header_info = re.match(r'^(\S+)\s+(\S+):(\d+)-(\d+)', header)
            if header_info:
                name = header_info.group(1)
                start = int(header_info.group(3))
                stop = int(header_info.group(4))
    
    def find_exons(self):
        'Adds exons to records information , exon is defined by strings of uppercase letters in fasta sequence'
        for record in self.records:
            sequence = record['sequence']
            exons = []
            i = 0
            while i < len(sequence):
                if sequence[i].isupper():
                    start = i
                    while i < len(sequence) and sequence[i].isupper():
                        i += 1
                    end = i
                    exons.append((start, end))
                else:
                    i += 1
            record['exons'] = exons
        return self.records


class Motif:
    def __init__(self, filename):
        self.filename = filename
        self.motifs = []
        self.patterns = []
    
    def get_motifs(self):
        'Fetches motifs from the motif file, storing in a list, implies that the format is one motif per line in file. '
        with open(self.filename) as f:
            for line in f:
                motif = self.motifs.append(line.strip())
            print(self.motifs)

    def convert_to_regex(self):
        'Converts a list of motifs to their corresponding regex based on nucleic acid notation'
        for motif in self.motifs:
            pattern = ''
            for char in motif:
                #Rules for nucleic acid notation.
                if char in ['a', 'A']:
                    pattern += '[Aa]'
                elif char in ['t', 'T']:
                    pattern += '[Tt]'
                elif char in ['c', 'C']:
                    pattern += '[Cc]'
                elif char in ['g', 'G']:
                    pattern += '[Gg]'
                elif char in ['y', 'Y']:
                    pattern += '[CcTt]'
                elif char in ['r', 'R']:
                    pattern += '[AaGg]'
                elif char in ['w', 'W']:
                    pattern += '[AaTt]'
                elif char in ['s', 'S']:
                    pattern += '[CcGg]'
                elif char in ['k', 'K']:
                    pattern += '[GgTt]'
                elif char in ['m', 'M']:
                    pattern += '[AaCc]'
                elif char in ['d', 'D']:
                    pattern += '[AaGgTt]'
                elif char in ['v', 'V']:
                    pattern += '[AaCcGg]'
                elif char in ['h', 'H']:
                    pattern += '[AaCcTt]'
                elif char in ['n', 'N']:
                    pattern += '[AaCcGgTt]'
                else:
                    pattern += char
            self.patterns.append(pattern)
        print(self.patterns)

    def search_motifs(self, records):
        'Returns a list of tuples for each motif match per record'
        results = []
        for record in records:
            for i, pattern in enumerate(self.patterns):
                matches = re.finditer(pattern, record['sequence'], re.IGNORECASE)
                for match in matches:
                    start = match.start() + 1  # base 1 conversion
                    stop = match.end()
                    #This is the format of the list. Record = header , index is order, start and stop are positions.
                    result = {'record': record, 'pattern_index': i, 'start': start, 'stop': stop}
                    results.append(result)
        return results
        

class Draw:
    #Define size and spacing for elements, color_map is the colors per motif
    def __init__(self, records, color_map=None):
        self.records = records
        self.width = 800
        self.height = 600
        self.font_size = 16
        self.line_spacing = 20
        self.rect_height = 5
        self.color_map = color_map or []

    def draw_motif_matches(self, results, filename):
        'Takes on positional data for motif and exon matches per record, and draws them to a .png with reference to record lines'
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        ctx = cairo.Context(surface)

        # Sets background to white
        ctx.set_source_rgb(1, 1, 1)
        ctx.paint()

        # Draws record lines
        ctx.set_line_width(2)
        ctx.set_source_rgb(0, 0, 0)
        for i, record in enumerate(self.records):
            y = ((i + 1) * self.height) / (len(self.records) + 2)
            ctx.move_to(0, y)
            ctx.line_to(self.width, y)
            ctx.stroke()

            # Adds header label per record line
            ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            ctx.set_font_size(self.font_size)
            text = record['header']
            x_bearing, y_bearing, width, height = ctx.text_extents(text)[:4]
            ctx.move_to((self.width - width) / 2, y - self.font_size - 2)
            ctx.show_text(text)

            # Draw motifs based on color_map
            ctx.save()
            for result in results:
                if result['record'] == record:
                    pattern_index = result['pattern_index']
                    #Currently defaults extra motifs to black
                    color = self.color_map[pattern_index] if pattern_index < len(self.color_map) else (0, 0, 0)
                    ctx.set_source_rgb(*color)
                    x = (result['start'] - 1) * self.width / len(record['sequence'])
                    y = ((i + 1) * self.height) / (len(self.records) + 2) + self.rect_height 
                    width = (result['stop'] - result['start'] + 1) * self.width / len(record['sequence'])
                    ctx.rectangle(x, y, width, self.rect_height)
                    ctx.fill()
            ctx.restore()

            # Draw Exons as Orange
            ctx.save()
            ctx.set_source_rgb(1, 0.5, 0)  # Set the color to orange
            for exon_start, exon_stop in record['exons']:
                exon_width = (exon_stop - exon_start) * self.width / len(record['sequence'])
                exon_x = exon_start * self.width / len(record['sequence'])
                exon_y = y - self.line_spacing - self.font_size - self.rect_height + 33
                ctx.rectangle(exon_x, exon_y, exon_width, self.rect_height)
                ctx.fill()

            ctx.restore()

                

        # Add legend for color map
        legend_x = self.width * 0.8
        legend_y = self.height * 0.8
        legend_width = self.width * 0.15
        legend_height = self.line_spacing * (len(self.color_map) + 1)
        ctx.rectangle(legend_x, legend_y, legend_width, legend_height)
        ctx.set_source_rgb(0, 0, 0)
        ctx.stroke()

        ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(self.font_size)
        for i, color in enumerate(self.color_map):
            y = legend_y + (i + 0.5) * self.line_spacing
            ctx.set_source_rgb(*color)
            ctx.rectangle(legend_x + 5, y - self.rect_height / 2, self.rect_height, self.rect_height)
            ctx.fill()
            ctx.move_to(legend_x + self.rect_height + 10, y)
            #Motif 1...2...3.. etc. refers to their order in the initial motif list.
            ctx.show_text(f"Motif {i+1}")

        # Add exon to the legend
        y = legend_y + (len(self.color_map) + 0.5) * self.line_spacing
        ctx.set_source_rgb(1, 0.5, 0)
        ctx.rectangle(legend_x + 5, y - self.rect_height / 2, self.rect_height, self.rect_height)
        ctx.fill()
        ctx.move_to(legend_x + self.rect_height + 10, y)
        ctx.show_text("Exon")

        surface.write_to_png(filename)




parser = FastaParse(fasta_file)
records = parser.parse_fasta()
parser.header_info()
parser.find_exons()
motif = Motif(motif_file)
motif.get_motifs()
motif.convert_to_regex()

matches = motif.search_motifs(records)

#To add more colors, simply add more pixel values to color_map below
color_map = [(1, 0, 0), (0, 1, 0), (0, 0, 1),(1,0.5,1)]
draw = Draw(records,color_map)

#Outputs drawing based on initial fasta file name
draw.draw_motif_matches(matches, f'{fasta_file}.png')

