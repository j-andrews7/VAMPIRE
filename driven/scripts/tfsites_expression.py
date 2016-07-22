#!/usr/bin/env python


import sys
import argparse
parser = argparse.ArgumentParser(usage=__doc__)
from math import ceil
from math import log10

####-CLASSES-####
# Not strictly necessary, but may be useful if other options are added
class Options_list:
    def __init__(self):
        #Should lines in the vcf output file be excluded 
        # if they don't have expression data?
        # -fe tag sets this to True
        self.filter_vcf = False

####-FUNCTIONS-####
def get_next_var( opened_file ):
    """
    Reads in the next line of the vcf and returns the next variant's info
    
    Args: opened_file = an already open input .vcf file
    
    Returns: None or a line info tuple with the following information
        (str list) motif names
        (str list) motif variant sequence scores
        (str list) motif reference sequence scores
        (boolean list) does the motif match overlap with a ChIP peak
        (str list) other fields
        (str) original line
        """
    
    line = opened_file.readline()
    
    line = line.strip()
    
    #input file is empty
    if line == "":
        return None
    
    #Find list of motifs that match
    line_list = line.split("\t")
    info_fields = line_list[7].split(";")
    #Holds other fields (other than motif fields which may be modified)
    other_info = []
    (motifns, motifvs, motifrs, motifcs) = ([],[],[],[])
    
    for field in info_fields:
        if field.startswith("MOTIFN="):
            #Cut off 'MOTIFN=' then convert to an array
            motifns = field[7:].split(',')
        elif field.startswith("MOTIFV="):
            motifvs = field[7:].split(',')
        elif field.startswith("MOTIFR="):
            motifrs = field[7:].split(',')
        #MOTIFC info field may be absent even if others are present
        elif field.startswith("MOTIFC="):
            motifcs = field[7:].split(',')
        else:
            other_info.append(field)
        
            

    return (motifns, motifvs, motifrs, motifcs, other_info, line)

def get_genes( opened_file, cell_line, threshold ):
    """
    Reads in and parses the .bed expression file
    
    Args: 
        opened_file = an already open input .bed file
        cell_line = (string) name of the cell line in question
        threshold = (float) expression threshold to count the gene as expressed
        
    Returns: gene_dict = (str -> float) a dictionary of:
        gene names -> expression levels
    """
    
    #find which column is the cell_line (contains relevant expression data)
    headers = opened_file.readline().strip().split('\t')
    cell_line_idx = -1
    for idx in range(len(headers)):
        if headers[idx].upper() == cell_line.upper():
            cell_line_idx = idx
            break
            
    if cell_line_idx == -1:
        print("Err: cell line <"+cell_line+"> not found in input file header.")
        return []
    
    gene_names = []
    exp_values = []
    
    #Read in expression levels for each gene
    for line in opened_file:  
        line_list = line.strip().split('\t')
        #Only add the gene to the list if its expression is above the threshold
        gene_name = line_list[3]
        exp_level = float(line_list[cell_line_idx])
        if exp_level >= threshold:
            #Add gene names and expression values to respective lists
            gene_names.append(std_gene_name(gene_name))
            exp_values.append(exp_level)
        
    return dict(zip(gene_names, exp_values))

def update_vcf(line_tup, output_f, options):
    """
    Updates the output file with the output in the correct variant call format
    
    Args:
        None or a line info tuple with the following information
            (str list) motif names
            (str list) motif variant sequence scores
            (str list) motif reference sequence scores
            (boolean list) does the motif match overlap with a ChIP peak
            (str list) other fields
            (str) original line
        output_f = output vcf
        options = options list
    
    Returns: Nothing (updates output_f instead of returning)
    """
    
    (motifns, motifvs, motifrs, motifcs, motifes, oth_info, line) = line_tup
    
    #First 8 columns should always be:
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    line = line.strip()
    columns = line.split('\t')
    
    #ID=MOTIFN,Type=String,Description="Matched motif names"
    #ID=MOTIFV,Type=Float,Description="Motif variant match scores"
    #ID=MOTIFR,Type=Float,Description="Motif reference match scores"
    #ID=MOTIFC,Type=Character,Description="Motif validated by ChIP (Y/N)"
    names = ""
    varscores = ""
    refscores = ""
    chips = ""
    explevels = ""
   
    
    for idx in range(len(motifns)):
        if names != "":
            names += ","
            varscores += ","
            refscores += ","
            chips += ","
            explevels += ","
        names += motifns[idx]
        varscores += motifvs[idx]
        refscores += motifrs[idx]
        explevels += sf_str(motifes[idx],4)
        if len(motifcs) != 0:
            chips += motifcs[idx]
    
    #If there are no matches, print the line unchanged or filter it out (return
    #without printing)
    if len(motifns) == 0 and options.filter_vcf == True:
        return
    
    outline = ""
    idx = 0
    
    for col in columns:
        if outline != "":
            outline += "\t"
        if idx == 7:
            if len(motifns) != 0:
                outline += "MOTIFN="+names+";MOTIFV="+varscores
                outline += ";MOTIFR="+refscores+";MOTIFE="+explevels
                if len(motifcs) != 0:
                    outline += ";MOTIFC="+chips
                for field in oth_info:
                    outline += ";" + field
            else:
                info = ""
                for field in oth_info:
                    if info != "":
                        info += ","
                    info += field
                outline += info
        else:
            outline += col
        idx += 1
    
    print (outline, file=output_f)
    
    return

def sf_str( x, n ):
    """ Rounds x to a certain number of significant figures.
        Returns the output as a string
        Ex:
        round(0.01234567, 3) -> 0.0123
        round(234.5678, 4) -> 234.6
        round(1234.56, 2) -> 1200 
        
    Args:
        x = float to be rounded
        n = number of sig figs
    
    Returns: a string """
    return str(round(x, int(n-ceil(log10(abs(x)))))) 

def std_gene_name(gene_name):
    """
    Returns standard format for gene name
        This function can be updated later to convert between standard names
    """
    return gene_name.upper()


####-PARSER-####  
# Create arguments and options
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-e", "--expression", dest = "exp_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-c" "--cell_line", dest = "cell_line", required=True)
parser.add_argument("-th", "--threshold", dest = "threshold", 
    required = False, default = 5)
parser.add_argument("-fe", "--filter_e", action="count", required = False)

args = parser.parse_args()

# Easier to use argument variables
inp_file = args.input_file
exp_file = args.exp_file
out_file = args.output_file
cell_line = args.cell_line
th = float(args.threshold)

#Options list. Easier to pass in functions or use in code updates.
options = Options_list()
options.filter_vcf = (args.filter_e != None)

####-MAIN-#### 
# Read in expression .bed file
gene_dict = []
with open(exp_file) as opened_ef:
    gene_dict = get_genes( opened_ef, cell_line, th )
#debug print(gene_dict)
if len(gene_dict) == 0:
    #break or something
    print("Error, no genes above threshold found")
  
# Open output file.
output_f = open(out_file,"w")


with open(inp_file) as vcf:

    line = vcf.readline()
    
    info_needed = True
    info = "##INFO=<ID=MOTIFE,Number=.,Type=Float,Description="
    info += "\"Motif expression level\">"
    
    # Skip info lines
    while line.startswith("##"):
    # Print new info lines at the top of the ##INFO section
        if info_needed and line.startswith("##INFO"):
            print(info, file=output_f)
            info_needed = False
        print(line, file=output_f, end="")
        line = vcf.readline()
    
    # First non-## line is a header
    print(line, file=output_f, end="")
    
    current_var = get_next_var( vcf )
    
    while current_var != None:
        
        (motifns, motifvs, motifrs, motifcs, oi, line) = current_var
        
        # Filter motifs by expression
        f_motifns = []
        f_motifvs = []
        f_motifrs = []
        f_motifcs = []
        # Expression level of motifs
        f_motifes = []
        
        for idx in range(len(motifns)):
            motif_name = motifns[idx]
            #Genes are only in the dictionary if their expression is above
            # the given threshold. Add them to the filtered lists if they are
            # expressed.
            #Standard gene name used (works if cases are different for now)
            std_name = std_gene_name(motif_name)
            if std_name in gene_dict:
                # Add expression level
                f_motifes.append(gene_dict[std_name])
                # Add name and other data
                f_motifns.append(motifns[idx])
                f_motifvs.append(motifvs[idx])
                f_motifrs.append(motifrs[idx])
                if (len(motifcs) != 0):
                    f_motifcs.append(motifcs[idx])
                    
        # Output filtered motifs
        tup = (f_motifns, f_motifvs, f_motifrs, f_motifcs, f_motifes, oi, line)
        update_vcf(tup, output_f, options)
        
        current_var = get_next_var( vcf )
    
output_f.close()
