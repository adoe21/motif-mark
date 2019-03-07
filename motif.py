import argparse
import cairo
import re
import random

#Collects all of the necessary inputs
def get_arguments():
	parser = argparse.ArgumentParser(description="Program to mark motifs on given sequences", add_help="False")
	parser.add_argument("-f", "--FASTA", help="FASTA file to be checked for motifs", required=True, type=str)
	parser.add_argument("-m", "--motif", help="File containing known motifs", required=True, type=str)
	return parser.parse_args()

#intiates objects for upcoming processes	
args = get_arguments()
stop = "start"
motif_dict = {}
og_motif = {}

#Creates a dictionary of motifs from motif input files
with open(args.motif, "r") as motifs:
	while stop != "":
		motif = motifs.readline()
		motif = motif.strip()
		if motif != "":
			motif_len = len(motif)
			#The following modifications and replace statements turn the motif into a regex object so Python can find the motifs in a sequence
			add_motif = '('+motif.upper()+')'
			Y_check = add_motif.replace('Y','[CT]')
			R_check = Y_check.replace('R','[AG]')
			S_check = R_check.replace('S','[GC]')
			W_check = S_check.replace('W','[AT]')
			K_check = W_check.replace('K','[GT]')
			M_check = K_check.replace('M','[AC]')
			B_check = M_check.replace('B','[CGT]')
			D_check = B_check.replace('D','[AGT]')
			H_check = D_check.replace('H','[ACT]')
			V_check = H_check.replace('V','[ACG]')
			N_check = V_check.replace('N','[ACGT]')
			U_check = N_check.replace('U','T')
			gap_check = U_check.replace('.', ' ')
			gap2_check = gap_check.replace('-',' ')
			motif_dict[gap2_check] = motif_len     #Creates a dictionary with the regex object as the key and the length of the original motif as the key
			og_motif[gap2_check] = add_motif       #Creates a dictionary with the regex object as the key and the original motif as the value
		stop = motif

#Creates a dictionary with the regex object as the key and rgb(red, green, blue) values as a key for future plotting of the motifs
motif_color = {}
for elem in motif_dict:
	motif_color[elem] = (random.uniform(0,1),random.uniform(0,1),random.uniform(0,1))

#Initiates more objects for upcoming processes	
stop="start"
sequence_dict = {}
longest_seq = 0
seq_motif = {}
seq_exon = {}

#Creates a dictionary from the input FASTA file which has the header of each sequence as the key and the sequence itself as the value
with open(args.FASTA, "r") as fasta:
	while stop != "":
		item = fasta.readline()
		stop = item.strip()
		if item.startswith(">")==True:
			header = item.strip()
			sequence_dict[header] = ""
		if item.startswith(">")==False:
			sequence = item.strip()
			sequence_dict[header] += sequence
			
				
			
for header in sequence_dict:            #Loops through all of the sequences 
	motif_mark = {}
	sequence = sequence_dict[header]
	if len(sequence) > longest_seq:
		longest_seq = len(sequence)
	if re.search(r'([A-Z]+)',str(sequence)) is None:    #Searches for exons and gives instructions for what to do if no exon is present
		for elem in motif_dict:
			motif_mark[elem] = []
			for match in re.finditer(elem,str(sequence)):   #Looks for motifs and saves their start position in a dictionary as the value with the motif as the key
				motif_mark[elem].append(match.start())
	if re.search(r'([A-Z]+)',str(sequence)) is not None:     #Searches for exons and gives instructions for what to do if an exon is present
		exon_start = re.search(r'([A-Z]+)',str(sequence)).span(1)[0]    #Records start and end position of exon
		exon_end = re.search(r'([A-Z]+)',str(sequence)).span(1)[1]
		for elem in motif_dict:
			motif_mark[elem] = []
			for match in re.finditer(elem,str(sequence).upper()):      #Again records motifs
				motif_mark[elem].append(match.start())
	seq_motif[header] = motif_mark
	seq_exon[header] = (exon_start, exon_end)              
#The two lines above record the motifs and exons for this specific sequence in dictionaries where the keys are the header of the sequence and the values are the motif start
#position dictionary and exon start and end positions.
	
surface = cairo.SVGSurface("plot.svg", longest_seq + 200, 200*len(sequence_dict))    #establishes the surface for plotting
context = cairo.Context(surface)
context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
start_height = 25
for seq in sequence_dict:
	context.set_source_rgb(0,0,0)
	context.set_line_width(1)
	context.move_to(50,start_height)
	context.set_font_size(20)
	context.show_text(seq)
	start_height += 50
	context.move_to(50,start_height)
	context.line_to(len(sequence_dict[seq])+50,start_height)
	context.stroke()
	context.move_to(seq_exon[seq][0]+50,start_height)
	context.set_line_width(20)
	context.line_to(seq_exon[seq][1]+50,start_height)
	context.stroke()
	#the above part of this loop draws the sequence with the intron symbolized by a thin line and the exon symbolized by a thick line
	for elem in seq_motif[seq]:
		context.set_source_rgb(motif_color[elem][0],motif_color[elem][1], motif_color[elem][2])
		for start_pos in seq_motif[seq][elem]:
			context.move_to(start_pos+50,start_height)
			context.set_line_width(20)
			context.line_to(start_pos+motif_dict[elem]+50,start_height)
			context.stroke()
		#The above part of the loop marks the motif on the sequence symbolized by a thick line with a color corresponding to which motif is being marked
	start_height += 100
color_start = 70
context.set_source_rgb(0,0,0)
context.move_to(longest_seq+50,color_start)
context.set_font_size(15)
context.show_text("Motif Color Legend:")
for elem in motif_color:
	color_start += 20
	context.set_source_rgb(motif_color[elem][0],motif_color[elem][1], motif_color[elem][2])
	context.move_to(longest_seq+50,color_start)
	context.show_text(og_motif[elem])
#The above creates a legend for the image with the original motif being listed in the color that it is marked with on the sequence