#!/usr/bin/env python 

import argparse
import math
import cairo
import re
import seaborn as sns

# Get args
parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",help="Input fasta file, ouput image will have same name and location.",required=True)
parser.add_argument("-m","--motif_fn",help="Input motifs file, up to 7 motifs.",required=True)
parser.add_argument("-l", "--lightmode", help="switch output to lightmode figure", action="store_true")
parser.add_argument("-s", "--size", nargs="?",type=int,const=1,default=2,help="scale factor to lengthen each nucleotide, default is 1")
args = parser.parse_args()

#########    Global Variables    #########

if not args.lightmode:
    # Use colorblind friendly palette from the seaborn package
    color_pal = [(0.00392156862745098, 0.45098039215686275, 0.6980392156862745), (0.8705882352941177, 0.5607843137254902, 0.0196078431372549), (0.00784313725490196, 0.6196078431372549, 0.45098039215686275), (0.8352941176470589, 0.3686274509803922, 0.0), (0.8, 0.47058823529411764, 0.7372549019607844), (0.792156862745098, 0.5686274509803921, 0.3803921568627451), (0.984313725490196, 0.6862745098039216, 0.8941176470588236), (0.9254901960784314, 0.8823529411764706, 0.2), (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]
    # Sequence and header color
    line_col = (.95,.95,.95)
else:
    # light mode colors
    color_pal = [(0.00392156862745098, 0.45098039215686275, 0.6980392156862745), (0.8705882352941177, 0.5607843137254902, 0.0196078431372549), (0.00784313725490196, 0.6196078431372549, 0.45098039215686275), (0.8352941176470589, 0.3686274509803922, 0.0), (0.8, 0.47058823529411764, 0.7372549019607844), (0.792156862745098, 0.5686274509803921, 0.3803921568627451), (0.984313725490196, 0.6862745098039216, 0.8941176470588236), (0.9254901960784314, 0.8823529411764706, 0.2), (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]
    line_col = (.4,.4,.4)

def draw_legend(c_key,x_pos,y_pos):
    x = x_pos+10
    y = y_pos+50
    context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(18)
    context.set_source_rgb(line_col[0],line_col[1],line_col[2]) 
    context.move_to(x,y-30)
    context.show_text("Legend:")
    for i in c_key:

        context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(18)
        context.set_source_rgb(c_key[i][0],c_key[i][1],c_key[i][2])
        context.move_to(x,y)
        context.show_text(i)
        x+=100
# width for exons
exon_width = 50

# introns line width
intron_width = 2

pad_header = 50
left_margin = 10

##### Define file parsing functions ######
   
def parse_fasta(file):
    '''Parses fasta file and returns a list of sequences and a list of fasta headers'''
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



#########    Define Class    ##########

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
        '''Draws the sequence'''
        # left margin
        x=self.origin[0] +5
        
        # space between seq line and fasta header
        y=self.origin[1]

        # arg to let user pick introns line width
        intron_width = 2
        
        # draw seq fasta header
        context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(16)
        context.set_source_rgb(line_col[0],line_col[1],line_col[2])
        context.move_to(x,y)
        context.show_text(self.header)

        for i in self.ex_in:            
            if i.group().isupper(): # check if exon or intron
                
                context.set_line_width(exon_width+args.size+4)
            else:
                
                context.set_line_width(intron_width)
            # Draw seq, with adjusted line width
            context.move_to(args.size*i.span()[0]+x,y+pad_header) 
            context.line_to(args.size*i.span()[1]+x,y+pad_header)
            context.set_source_rgb(line_col[0],line_col[1],line_col[2]) 
            context.set_dash([1,0])
            context.stroke()
            

        
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
            # Find the motifs in the sequence, returns zero length matches
            match_res = re.finditer("(?=("+"".join(regex_list)+"))",self.seq.upper())
            for i in match_res:
                match_dict[m].append(i.span()) 
        
        jitter =-1
        jit_dict = {} # Used to offset the y axis of each motif
        for i in match_dict:
            jitter+=.5
            jit_dict[i]=jitter # add .5 in increments and assign to each motif
            
            for k in match_dict[i]:

                # height for start point
                rect_start_y=(y+pad_header)-exon_width/2 +jit_dict[i]
                # draw a rectangle:
                # Top left
                context.move_to(args.size*k[0]+x,rect_start_y)
                # Top right
                len_m = (args.size*k[0]+x+len(i)*args.size)-(args.size/2) # account for outline width
                context.line_to(len_m,rect_start_y) 
                # Bottom right
                context.line_to(len_m,rect_start_y+exon_width)
                # Bottom left
                context.line_to(args.size*k[0]+x,rect_start_y+exon_width)
                # Top left
                context.line_to(args.size*k[0]+x,rect_start_y)
                context.close_path()
                # fill rectangle
                context.set_source_rgba(col_key[i][0],col_key[i][1],col_key[i][2],.25) 
                context.fill_preserve()
                # outline
                context.set_source_rgba(col_key[i][0],col_key[i][1],col_key[i][2],.85)  
                context.set_line_width(args.size)
                context.stroke()   

#########    Parse inputs     #########

# Parse input file
file_lists = parse_fasta(args.filename)
headers = file_lists[1]
seqs = file_lists[0]

# Parse motif file
motifs = parse_motifs(args.motif_fn)

#########    Create Objects     ##########

# dictionary to track colors for each motif
col_key ={} 
col_count = 0
for i in motifs:
    col_key[i]=color_pal[col_count]
    col_count+=1

ob_dict = {}
start_y = 20 # starting y pos

for i in range(len(headers)):

    ob_dict[i]=motif_mark((left_margin,start_y),seqs[i],motifs,headers[i])
    i+=1
    start_y+=125 # increment by 100 to move the y origin for each figure
    # *could add an option here to incerase the spacing between plots* 

#########    Draw Background     ##########


# width = max width of each figure + padding for right margin 
# height = number of motif_marks * width of each + padding for legend
WIDTH, HEIGHT = len(max(seqs))*args.size+250*args.size, len(seqs)*120+150
# create the coordinates to display the graphic
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, WIDTH, HEIGHT)

# create the coordinates to draw on
context = cairo.Context(surface)


# fill background
if not args.lightmode:
    context.set_source_rgb(0.1, 0.1, 0.1)
else:
    context.set_source_rgb(1,1,1)
context.rectangle(0, 0, WIDTH, HEIGHT)
context.fill()

####### Draw figures using the objects ######
for i in ob_dict:
    ob_dict[i].draw_seq()
    ob_dict[i].draw_motifs()

# Draw legend
draw_legend(col_key,left_margin,start_y)


#write out drawing
surface.write_to_png (args.filename.split(".")[0]+".png") # Output to PNG