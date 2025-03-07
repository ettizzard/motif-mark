#!/usr/bin/env python

# Importing necessary packages
from itertools import product
import re
import cairo
import argparse 


# Defining Functions

# Initializing argparse
def get_args():
    parser = argparse.ArgumentParser(description="Script which generates a .png plot of motifs on genes given a motif and fasta file.")
    parser.add_argument("-f", "--fastafile", help = "input fasta file with gene names in header", type = str)
    parser.add_argument("-m", "--motiffile", help = "input motif text file with one motif per line", type = str)


    return parser.parse_args()

args = get_args()


# In the event that the input fasta file is not properly formatted, the sequence lines will need to be concatenated to a single line.
def generate_list_of_fasta_lines(fasta_file_input) -> list:
    '''This function is basically just a oneline fasta function, but instead appends each line to a list rather than to a new file.'''
    first_line = True
    list_of_fasta_lines = []
    
    with open(fasta_file_input,"r") as fasta_file:
        sequence_line = ""
       
        for line in fasta_file:
            line = line.strip()
           
            if line.startswith(">"):
                if first_line == True:
                    list_of_fasta_lines.append(line)
                    first_line = False
                
                else:
                    list_of_fasta_lines.append(sequence_line)
                    list_of_fasta_lines.append(line)
                    
                    sequence_line = ""
            else:
                sequence_line += line
    
    list_of_fasta_lines.append(sequence_line)
    
    return list_of_fasta_lines



def generate_all_gene_sequences_from_motif(input_motif):
   '''Takes in an ambiguous motifs and outputs all possible gene sequences which follow that motif. If motif in unambiguous, returns just the unmodified sequence.'''
   degenerate_bases = {"R": ["A", "G"], 
                    "Y": ["C", "T"], 
                    "S": ["C", "G"], 
                    "W": ["A", "T"], 
                    "K": ["G", "T"], 
                    "M": ["A", "C"], 
                    "B": ["C", "G", "T"], 
                    "D": ["A", "G", "T"], 
                    "H": ["A", "C", "T"], 
                    "V": ["A", "C", "G"], 
                    "N": ["A", "C", "G", "T"]}
   
   # Adding each character to list individually. If character is a degenerate base, all of its possible bases are added
   list_of_str_characters = [(character,) if character not in degenerate_bases else degenerate_bases[character] for character in input_motif]
   # Using product function, making the cartesian products of our characters to produce all possible gene sequences from the motif
   return (''.join(o) for o in product(*list_of_str_characters))


def generate_set_of_all_possible_motifs(motif) -> set:
    '''Takes in the possible sequences from the above function, enforces upper capitalization, and adds them to a set.'''
    motif = motif.upper()
    
    set_of_all_possible_motifs = set(generate_all_gene_sequences_from_motif(motif))    

    return set_of_all_possible_motifs


def find_motif_positions(query_motif_set:set, gene_sequence:str) -> set:
    '''Takes a set of all possible sequences represented by a motif as well as a gene sequence as input. The positions of the motif sequences within the gene are then output in a new set of tuples.'''
    set_of_motif_positions = set()
    
    # Search for all instances of motifs
    for motif in query_motif_set:
        # Utilize a lookahead expression to find overlapping motifs as well
        matches = re.finditer(f'(?=({motif}))', gene_sequence)
    
        # For every motif found, add its start and end position to the set
        for match in matches:
            set_of_motif_positions.add((match.start(1), match.end(1)))
    
    # Resulting set looks like: {(0, 3), (7, 10)...}
    return set_of_motif_positions




# For plotting, we need to determine the type of feature, feature sequence, length, and location of all gene features in each fasta gene.
def split_gene_feature_regions(input_gene_string: str) -> list:
    '''Takes an input gene string which utilizes capitalization to delineate exonic and intronic regions and outputs a list of the gene's feature type, sequence, length, and starting position.'''
    # Boolean checking if gene sequence starts in an exonic region - True:Exon, False:Intron
    current_feature_is_exonic = input_gene_string[0].isupper()
    
    list_of_feature_locations = []
    current_feature = ""
    position_counter = 0
    
    # For every nucleotide in the gene string, check for exonic/intronic identity based on capitalization.
    for nucleotide in input_gene_string:
        
        # If current nucleotide is same type as current feature, concatenate nucleotide to current feature string, increment position, and calculate feature starting position.
        if nucleotide.isupper() == current_feature_is_exonic:
            current_feature += nucleotide
            position_counter += 1
            feature_starting_position = position_counter - len(current_feature)
        
        # Otherwise, a new feature type has been encountered, and append the finished feature to the list of all feature data. Then, reset current feature
        # string to current nucleotide, increment position, and set the new feature type checker to the opposite boolean value.
        else:
            if nucleotide.isupper() == False:
                list_of_feature_locations.append(("Exon", current_feature, len(current_feature), feature_starting_position))
            
            else:
                list_of_feature_locations.append(("Intron", current_feature, len(current_feature), feature_starting_position))
                
            current_feature = nucleotide
            position_counter += 1
            current_feature_is_exonic = not current_feature_is_exonic
    
    # Final append statements for last feature that needs to be outside of original loop
    if nucleotide.isupper() == True:
        list_of_feature_locations.append(("Exon", current_feature, len(current_feature), feature_starting_position))
    
    else:
        list_of_feature_locations.append(("Intron", current_feature, len(current_feature), feature_starting_position))
    
    
    return list_of_feature_locations


# Setting up files

# Declaring variable to contain list of fasta file lines for accession
input_fasta = args.fastafile
fasta_file = generate_list_of_fasta_lines(input_fasta)

# Opening motif file as variable for later access
input_motif = args.motiffile
motif_file = open(input_motif,'r')

# Setting name for output file
generated_plot_filename = f"{input_fasta.split('.f')[0]}.png"


# Defining classes

class Gene:
    def __init__(self, header_line, sequence_line):
        # Removing all characters but gene name from header line
        self.name = (header_line.strip().split(" ")[0]).split(">")[1]
        
        # Gene sequence is just fasta file sequence line
        self.sequence = sequence_line.strip()
        
        # Generating list of gene feature regions
        self.list_of_introns_and_exons = split_gene_feature_regions(self.sequence)
        
        # When searching gene sequence for motifs, everything needs to be consistent in capitalization
        self.sequence_uppercase = self.sequence.upper()
        
        # We need length of the total sequence for plotting scale
        self.sequence_length = len(self.sequence)





# Setting up color palette for plotting
#                 orange          green         blue             purple         pink    
color_palette = [(1, 0.412, 0), (0.149, 1, 0), (0, 0.859, 1), (0.643, 0, 1), (1, 0, 0.816)]



class Motif:
    # Setting up a counter to change color of motif object when a new motif is accessed
    motif_times_called = 0
    
    def __init__(self, motif_file_line):
        
        # Unmodified name is the motif sequence present in the input motif file. Since this could include motifs with 
        # inconsistent capitalization and could be RNA rather than DNA, another attribute will be needed for the normalized motifs.
        self.unmodified_name = motif_file_line.strip()

        # Converting any RNA to DNA
        if "u" in self.unmodified_name or "U" in self.unmodified_name:
            self.name = self.unmodified_name.replace('u','T')
            self.name = self.unmodified_name.replace('U','T')
        
        else:
            self.name = self.unmodified_name

        # Assigning color to each motif
        self.color = color_palette[Motif.motif_times_called]
        Motif.motif_times_called += 1
        
        self.length = len(self.name)
        
        # Assigning all possible sequences based on the motif
        self.sequences = generate_set_of_all_possible_motifs(self.name)
        


class MotifPerGene:
    motif_gene_name_set = set()
    
    def __init__(self, gene_object, motif):
        self.motif_name = motif.unmodified_name
        self.gene_name = gene.name
        self.motif_positions = find_motif_positions(motif.sequences, gene_object.sequence_uppercase)
        self.color = motif.color
        MotifPerGene.motif_gene_name_set.add(self.gene_name)




# Creating lists of class objects

# For every header line in fasta file, assign to variable. For every sequence line, assign to variable, then append list with
# gene object created from both variables
list_of_gene_objects = []

for line_num, line in enumerate(fasta_file):
    
    if line_num %2 == 0:
        header_line = line
    
    else:
        sequence_line = line
        list_of_gene_objects.append(Gene(header_line, sequence_line))



# For every line in motif file, create a motif object and append to list
list_of_motif_objects = []

for line in motif_file:
    list_of_motif_objects.append(Motif(line))        


# Now in the meta-listicon. Iterate through the gene list - for each gene, iterate through the motif list. For every motif in
# the motif object list, append the list of motif_per_gene objects with the MotifPerGene object created from that set of gene
# and motif.
list_of_motif_per_gene_objects = []

# Collecting length of all gene sequences in a list, then finding the maximum value for figure dimensions
list_of_gene_sequence_lengths = []

for gene in list_of_gene_objects:
    list_of_gene_sequence_lengths.append(len(gene.sequence))
    
    for motif in list_of_motif_objects:
        list_of_motif_per_gene_objects.append(MotifPerGene(gene, motif))



# Cairo code and loops

# Setting up dimensional variables
longest_gene_sequence_length = max(list_of_gene_sequence_lengths)

number_of_genes = len(MotifPerGene.motif_gene_name_set)

font_size = 40


# Setting up cairo surface and background

surface_width = round(longest_gene_sequence_length * 1.5)

# Making length of figure vary by number of genes provided 
surface_length = round(number_of_genes * 200)


surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, surface_width, surface_length)

context = cairo.Context(surface) 

# Setting white background
#context.set_source_rgb(1, 1, 1)
context.set_source_rgb(0.902, 0.878, 0.851)
context.rectangle(0, 0, surface_width, surface_length)
context.fill()


# Setting up black color and font size for printing gene names
context.set_source_rgb(0, 0, 0) 
context.set_font_size(font_size) 
context.select_font_face("Futura", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
context.set_antialias(cairo.ANTIALIAS_BEST)


# Plotting gene names
left_border = 50

# vertical_shift is the top border for the first gene, then is incremented for every gene after
top_border = 125

# Printing Figure title
context.move_to(left_border, 50)
context.show_text("Plot of all provided motifs on their respective genes")


gene_name_vertical_shift = -5

gene_label_heights = []
for gene in list_of_gene_objects:
    # Move brush to starting position within image border
    context.move_to(left_border, top_border + gene_name_vertical_shift)
    
    # Print gene name
    context.show_text(f'"{gene.name}" Gene') 
    
    # using text_extents function to fetch official height of gene name string in figure
    gene_name_height = (context.text_extents(f'"{gene.name}" Gene')[3])
    gene_label_heights.append(gene_name_height)
    
    gene_name_vertical_shift += 200
    context.stroke() 

#print(gene_label_heights)


# Plotting intronic and exonic lines for gene
intron_line_width = 8
exon_line_width = intron_line_width * 3

feature_vertical_shift = top_border + 35
gene_positions = {}
for gene in list_of_gene_objects:
    for feature in gene.list_of_introns_and_exons:
        
        if feature[0] == 'Intron':
            context.set_line_width(intron_line_width)
        else:
            context.set_line_width(exon_line_width)
        
        # feature[3] = feature starting position
        # Moving brush to starting position of motif
        context.move_to(left_border + feature[3], feature_vertical_shift)


        # feature[2] = length of feature
        # Drawing line from starting positon to ending position of motif
        context.line_to(left_border + feature[3] + feature[2], feature_vertical_shift)
        context.stroke()
    gene_positions[gene.name] = feature_vertical_shift
    feature_vertical_shift += 200



# Plotting motif lines against gene
for motif in list_of_motif_per_gene_objects:
    for gene in motif.motif_gene_name_set:
        if motif.gene_name == gene:
            
            for position in motif.motif_positions :

                context.set_source_rgba(motif.color[0], motif.color[1], motif.color[2], 0.5) 
                context.rectangle(left_border + position[0], gene_positions[gene] - 25, position[1] - position[0], 50)
                context.fill()



# Legend code
context.set_font_size(15)

legend_vertical_shift = 0

for motif in list_of_motif_objects:
    
    # Draw color blocks for legend
    context.set_line_width(15)
    context.set_source_rgba(motif.color[0], motif.color[1], motif.color[2], 0.5)
    context.move_to(surface_width - (left_border * 4), top_border + legend_vertical_shift)
    context.line_to(surface_width - (left_border * 3), top_border + legend_vertical_shift)
    context.stroke()
    
    # Print motif label per color block
    context.set_source_rgb(0, 0, 0) 
    context.move_to(surface_width - left_border - 100, top_border + legend_vertical_shift + 5)
    context.show_text(motif.unmodified_name)
    context.stroke()
    legend_vertical_shift += 20


# Intron legend section
# Print "intron" legend label 
context.set_source_rgb(0, 0, 0) 
context.move_to(surface_width - left_border - 100, top_border + legend_vertical_shift + 5)
context.show_text('intron')
context.stroke()
# Print intron line width
context.set_line_width(intron_line_width)
context.move_to(surface_width - (left_border * 4), top_border + legend_vertical_shift)       
context.line_to(surface_width - (left_border * 3), top_border + legend_vertical_shift)
context.stroke()

legend_vertical_shift += 20

# Exon legend section
# Print "exon" legend label
context.move_to(surface_width - left_border - 100, top_border + legend_vertical_shift + 5)
context.show_text('exon')
context.stroke()
# Draw exon line width
context.set_line_width(exon_line_width)
context.move_to(surface_width - (left_border * 4), top_border + legend_vertical_shift)       
context.line_to(surface_width - (left_border * 3), top_border + legend_vertical_shift)
context.stroke()

# Write to output png
surface.write_to_png(generated_plot_filename)

# Successful exit message
print(f"{generated_plot_filename} successfully generated :)")