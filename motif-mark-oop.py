#!/usr/bin/env python 

import argparse
import math
from unicodedata import numeric
import cairo
import re
import numpy as np
import seaborn as sns

# Get args
parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",help="Input fasta file",required=True)
parser.add_argument("-m","--motif_fn",help="Input motifs file",required=True)
parser.add_argument("-d", "--darkmode", help="ouputs dark mode figures", action="store_true")
parser.add_argument("-s", "--size", nargs="?",type=int,const=1,default=1,help="scale factor to lengthen each nucleotide")
args = parser.parse_args()


##### Define file parsing functions ######
   
def parse_fasta(file):
    with open(file, "r") as fh:
        
        header_list = []
        seq_list = []
        is_first = 0
        seq=""
        
        for line in fh:
            
            line = line.strip()
            
            if line.startswith(">"):
                if is_first == 0:
                    is_first+=1
                else:
                    seq_list.append(seq)
                    seq=""
                header_list.append(line)
            else:

                seq = seq+line
        seq_list.append(seq)
    return seq_list,header_list

def parse_motifs(file):
    with open(file, "r") as fh:
        motif_list = []
        for line in fh:
            line = line.strip()
            motif_list.append(line)
    return motif_list

# Dictionary for degenerate IUPAC bases
bases_dict = {"A":"A",
               "C":"C",
               "G":"G",
               "T":"T",
               "U":"U",
               "W":"[AT]",
               "S":"[CG]",
               "M":"[AC]",
               "K":"[GT]",
               "R":"[AG]",
               "Y":"[CT]",
               "B":"[CGT]",
               "D":"[AGT]",
               "H":"[ACT]",
               "V":"[ACG]",
               "N":"[ACGT]"}



#########    Define Classes     ##########

class motif_mark:
    def __init__(self, origin, seq, motif_list,header):
        '''motif_mark class for each individual figure'''
      ## Data ##
        self.origin = origin
        self.seq = seq
        self.motif_list = motif_list
        self.header = header

        self.seq_len = len(seq)

        # get exon and intron posititions and lengths in order
        self.ex_in = re.finditer("[A-Z]+|[a-z]+",self.seq)
        

    
  ## Methods ##
    def draw_seq(self):

        # left margin
        x=self.origin[0] +5
        
        # space between seq line and fasta header
        y=self.origin[1]

        # arg to let user pick introns line width
        intron_width = 2
        
        # draw seq fasta header
        context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(14)
        context.set_source_rgb(line_col[0],line_col[1],line_col[2])
        context.move_to(x,y)
        context.show_text(self.header)

        for i in self.ex_in:
            #print(i.span())
            #print(i.group())
            
            if i.group().isupper(): # check if exon or intron
                
                context.set_line_width(exon_width)
            else:
                
                context.set_line_width(intron_width)
                
            context.move_to(args.size*i.span()[0]+x,y+pad_header) 
            context.line_to(args.size*i.span()[1]+x,y+pad_header)
            context.set_source_rgb(line_col[0],line_col[1],line_col[2]) 
            context.stroke()
            
    def are_overlapping(r, s):
        return not(r[1] < s[0] or s[1] < r[0])
        
    def draw_motifs(self):
        # left margin
        x=self.origin[0] +5
        
        # space between seq line and fasta header
        y=self.origin[1]
        match_dict={}
        for m in self.motif_list:
            # key = motif, value = (start,end)
            match_dict[m] = []
            regex_list = [bases_dict[i] for i in m.upper()]
            
            match_res = re.finditer("".join(regex_list),self.seq.upper())
            for i in match_res:
                match_dict[m].append(i.span()) 
        
        
        for i in match_dict:
            for k in match_dict[i]:

                context.set_line_width(exon_width)
                context.move_to(args.size*k[0]+x,y+pad_header) 
                context.line_to(args.size*k[1]+x,y+pad_header)
                
                context.set_source_rgba(col_key[i][0],col_key[i][1],col_key[i][2],.7) 
                context.stroke()   


#########    Parse inputs     #########

# Parse input file
file_lists = parse_fasta(args.filename)
headers = file_lists[1]
seqs = file_lists[0]

# Parse motif file
motifs = parse_motifs(args.motif_fn)





#########    Global Variables    #########

if args.darkmode:
    # Get colorblind friendly palette
    color_pal = sns.color_palette("colorblind")
    # Remove Gray
    color_pal.remove((0.5803921568627451, 0.5803921568627451, 0.5803921568627451))
    line_col = (1,1,1)
else:
    # light mode colors
    color_pal = sns.color_palette("colorblind",desat=.5)
    # Remove Gray
    color_pal.remove((0.5803921568627451, 0.5803921568627451, 0.5803921568627451))
    line_col = (0,0,0)



col_key ={} #track colors for each motif
col_count = 0
for i in motifs:
    col_key[i]=color_pal[col_count]
    col_count+=1


# replace with arg to let user pick width for exons
exon_width = 50

# replace with arg to let user pick introns line width
intron_width = 2

pad_header = 50
left_margin = 10

#########    Create Objects     ##########

ob_dict = {}

start_y = 20


for i in range(len(headers)):
    ob_dict[i]=motif_mark((left_margin,start_y),seqs[i],motifs,headers[i])
    i+=1
    start_y+=100



#########    Draw Background     ##########


# width = max width of each figure 
# height = number of motif_marks * width of each + padding
WIDTH, HEIGHT = 1000+300*args.size, len(seqs)*100+200

# create the coordinates to display the graphic
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, WIDTH, HEIGHT)

# create the coordinates to draw on
context = cairo.Context(surface)


# fill background
if args.darkmode:
    context.set_source_rgb(0.23529411764705882, 0.23529411764705882, 0.23529411764705882)
else:
    context.set_source_rgb(1,1,1)
context.rectangle(0, 0, WIDTH, HEIGHT)
context.fill()

####### Draw using the objects ######


for i in ob_dict:
    ob_dict[i].draw_seq()
    ob_dict[i].draw_motifs()




#write out drawing
surface.write_to_png ("test.png") # Output to PNG