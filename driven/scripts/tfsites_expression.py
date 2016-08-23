#!/usr/bin/env python

"""
Last modified: 08/23/2016
Authors: Ethan Pfeifer & Jared Andrews

For a given motif annotated vcf file (already run through tfsites_checker.py)
remove all motif matches for TFs that are not expressed in the given sample.
    
Usage: tfsites_checker.py -i <input.vcf> -e <expression.bed>
 -o <output.txt> -c <sample_name>

Args:
    -i (required) <input.vcf> = Name of sorted variant file to process. 

    -o (required) <output.vcf> = Name of output file to be created.

    -e (required) <expression.bed> = An expression 'bed' file.

    -sn (optional) <sample_name> = Name of the sample (the name must be present
        in expression.bed).

    -th (optional) <5> = Motifs are expressed a match if they are expressed above
      a given threshold. 

    -fe (optional flag) = If used, variants that do not match any motif for 
      an expressed protein will not be included in the output (-o) file.
"""

#Note: currently this code expects an expression text file of the form:
# GENE_NAME <Sample names delineated by tabs>
#To change where gene name is expected, go to line 128
# gene_name = line_list[0]
# where [0] works for the above configuration and [3] works for a 'bed' file
# 'bed' file: chr	START	END	GENE_SYMB	<Sample names delineated by tabs>

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
        (list list) other motif info fields
            each other field has the info name as the first object
            ex: ["MOTIFV",["0.5", "0.24"]] -> MOTIFV=0.5;0.24
        (str list) other info fields
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
    (motifns, motif_other) = ([],[])
    
    for field in info_fields:
        if field.startswith("MOTIFN="):
            #Cut off 'MOTIFN=' then convert to an array
            motifns = field[7:].split(',')
        #Check for other motif fields such as:
        #(MOTIFN)   MOTIFV  MOTIFR  MOTIFC
        # MOTIFVH   MOTIFRH MOTIFVG MOTIFRG
        elif (field.startswith("MOTIFV=") or field.startswith("MOTIFR=") or
            field.startswith("MOTIFC=")):
            name = field[:7]
            values = field[7:].split(',')
            motif_other.append([name,values])
        elif (field.startswith("MOTIFVH=") or field.startswith("MOTIFRH=") or
            field.startswith("MOTIFVG=") or field.startswith("MOTIFRG=")):
            name = field[:8]
            values = field[8:].split(',')
            motif_other.append([name,values])
        else:
            other_info.append(field)
            

    return (motifns, motif_other, other_info, line)

def get_genes( opened_file, cell_line, threshold ):
    """
    Reads in and parses the .bed expression file
    
    Args: 
        opened_file = an already open input .bed file
        cell_line = (string) name of the cell line or sample in question
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
        print("Err: Sample <"+cell_line+"> not found in input file header.")
        return []
    
    gene_names = []
    exp_values = []
    
    #Read in expression levels for each gene
    for line in opened_file:  
        line_list = line.strip().split('\t')
        #Only add the gene to the list if its expression is above the threshold
        gene_name = line_list[0]
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
            (str list list) list of lists of other motif info
            (str list) other info fields
            (str) original line
        output_f = output vcf
        options = options list
    
    Returns: Nothing (updates output_f instead of returning)
    """
    
    (motifns, motifos, motifes, oth_info, line) = line_tup
    
    #First 8 columns should always be:
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    line = line.strip()
    columns = line.split('\t')
    
    #Strings for motif info fields
    str_m_os = ["" for idx in range(len(motifos))]
    
    #Stored as:
    # ["MOTIFC", ['Y','N','Y','Y']] -> MOTIFC=Y,N,Y,Y
    names = ""
    explevels = ""
    
    for idx in range(len(motifns)):
        if idx != 0:
            names += ","
            for idy in range(len(str_m_os)):
                str_m_os[idy] += ","
            explevels += ","
        names += motifns[idx]
        for idy in range(len(str_m_os)):
            str_m_os[idy] += motifos[idy][1][idx]
        explevels += sf_str(motifes[idx],4)
    
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
                outline += "MOTIFN="+names+";MOTIFE="+explevels
                for idy in range(len(motifos)):
                    outline += ";"+motifos[idy][0]+str_m_os[idy]
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
parser.add_argument("-sn" "--sample_name", dest = "sample_name",
    required=False, default = None)
parser.add_argument("-th", "--threshold", dest = "threshold", 
    required = False, default = 5)
parser.add_argument("-fe", "--filter_e", action="count", required = False)

args = parser.parse_args()

# Easier to use argument variables
inp_file = args.input_file
exp_file = args.exp_file
out_file = args.output_file
th = float(args.threshold)

#Options list. Easier to pass in functions or use in code updates.
options = Options_list()
options.filter_vcf = (args.filter_e != None)

#Get sample name
sample_name = args.sample_name
if sample_name == None:
    #Last part of input after slashes
    file_name = inp_file.split('\\')[-1].split('/')[-1]
    #First part of filename before period
    sample_name = file_name.split('.')[0]
    
####-MAIN-#### 
# Read in expression .bed file
gene_dict = []
with open(exp_file) as opened_ef:
    gene_dict = get_genes( opened_ef, sample_name, th )
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
    info += "\"Motif expression level z-score\">"
    
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
        
        (motifns, motifos, oi, line) = current_var
        
        # Filter motifs by expression
        f_motifns = []
        
        # Copy names into filtered motifs list
        f_motifos = []
        for idx in range(len(motifos)):
            #Stored as:
            # ["MOTIFC", ['Y','N','Y','Y']] -> MOTIFC=Y,N,Y,Y
            info_name = motifos[idx][0]
            f_motifos.append( [info_name, []] )
        
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
                for idy in range(len(motifos)):
                    f_motifos[idy][1].append(motifos[idy][1][idx])
                    
        # Output filtered motifs
        tup = (f_motifns, f_motifos, f_motifes, oi, line)
        update_vcf(tup, output_f, options)
        
        current_var = get_next_var( vcf )
    
output_f.close()
