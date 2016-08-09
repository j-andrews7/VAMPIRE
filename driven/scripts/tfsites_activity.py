#!/usr/bin/env python

import os
import sys
import argparse
#Standard deviation and mean for z-scores
from numpy import std
from numpy import mean
#Needed to round by significant figures
from math import ceil
from math import log10

####-CLASSES-####

#Vcf class is used to manage input / output of individual vcfs.
# Fields:
#   name = sample name for given vcf, used for printing
#   vcf_f = opened vcf file, only used internally
#   out_f = opened vcf file, only used internally (may be None)
#   var_pos = Position object
#   matches = enhancer activity matches
#   original_line = (string) original info line from vcf
# May throw FileNotFoundError
class Vcf:
    def __init__(self, sample_name, input_path, output_path):
        self.name = sample_name
        self.vcf_f = open(input_path)
        self.matches = []
        if output_path != None:
            self.out_f = open(output_path,"w")
        else:
            self.out_f = None
        
        #Initialize current variant (aka current line)
        line = self.vcf_f.readline()
        info_needed = True
        info = "##INFO=<ID=EZACT,Number=1,Type=Float,Description="
        info += "\"Enhancer z-score for activity "
        info += "(variant activity vs reference activity)\">"
        
        #Skip info lines
        while line.startswith("#"):
        #Print new info lines at the top of the ##INFO section
            if self.out_f != None:
                if info_needed and line.startswith("##INFO"):
                    print(info, file=self.out_f)
                    info_needed = False
                print(line, file=self.out_f, end="")
            line = self.vcf_f.readline()
        
        #Current variant
        self.parse_line(line)
        
    #Outputs old variant if output path is given
    def next_variant(self, options):
        if self.out_f != None:
            self.output_var(options)
        
        self.matches = []
        line = self.vcf_f.readline()
        
        #skip info lines-
        while line.startswith("#"):
            line = self.vcf_f.readline()

        self.parse_line(line)
    
    #Remember to check if self.out_f is present before calling output_var
    def output_var(self, options):
    
        #If the filter vcf option is on, don't print variants that don't
        # have an expression z-score (above the threshold if present)
        if options.filter_vcf and len(self.matches) == 0:
            return
        
        line_list = self.original_line.split("\t")
        
        #Initialize the line to be printed
        line = line_list[0]
        for idx in range(1,len(line_list)):
            #The INFO column
            if idx == 7 and len(self.matches) > 0:
                info = ""
                #Matches should have enhancer id, Position object, and zscore
                for match in self.matches:
                    (id, pos, zscore) = match
                    if info != "":
                        info += ";"
                    info += "EZACT="+sf_str(zscore,3)
                if line_list[idx] != '.':
                    info += ';'+line_list[idx]
                line += '\t'+info
            else:
                line += '\t'+line_list[idx]
        
        print(line, file=self.out_f)
    
    def parse_line(self, line):
        """
        Reads in the next line of the vcf and returns the next variant's information
        
        Args: line = (string) the vcf line to parse
        
        Returns: None or a tuple with the following information (in order):
            The position of the variant as a Position class object (see below)
            A list of the motif names that matched that variant
            A list of other motif fields for futher processing
            A list of other INFO fields
            The original line
        """
        line = line.strip()
        
        #input file is empty
        if line == "":
            self.var_pos = None
            #self.motif_names = None
            #self.motif_other_info = None
            #self.other_info_fields = None
            self.original_line = None
            return
        
        line_list = line.split("\t")
        
        #Find variant position
        chr = line_list[0]
        start_pos = int(line_list[1])
        ref_seq = line_list[3]
        end_pos = start_pos + len(ref_seq)
        var_pos = Position(chr, start_pos, end_pos)
        """
        #Find list of motifs that match
        info_fields = line_list[7].split(";")
        #Names of motifs that matched the variant (MOTIFN)
        motifns = []
        #Variant sequence motif match scores (MOTIFV)
        #motifvs = []
        #Reference sequence motif match scores (MOTIFR)
        #motifrs = []
        #Other motif info fields such as:
        # MOTIFC, MOTIFE
        motif_others = []
        #Holds other fields (other than motif fields which may be modified)
        others = []
        
        
        for field in info_fields:
            if field.startswith("MOTIFN="):
                #Cut off 'MOTIFN=' then convert to an array
                motifns = field[7:].split(',')
            #elif field.startswith("MOTIFV="):
            #    motifvs = field[7:].split(',')
            #elif field.startswith("MOTIFR="):
            #    motifrs = field[7:].split(',')
            #May eventually need to modify this to be more specific, or
            # tell the user not to use other INFO fields beginning with MOTIF
            elif field.startswith("MOTIF"):
                motif_info_list = field[7:].split(',')
                motif_others.append(motif_info_list)
            else:
                others.append(field)
        

        self.motif_names = motifns
        self.motif_other_info = motif_others
        self.other_info_fields = others
        """
        self.var_pos = var_pos
        self.original_line = line
        
        return
    
    def close(self):
        self.vcf_f.close()
        #self.out_f.close()
        
    def __str__(self):
        return self.name

# Fields:
#   chr = (string) Chromosome (chr1, chr2, etc)
#   start = start position
#   end = end position
class Position:
    def __init__(self, chr, start_pos, end_pos):
        self.chr = chr
        self.start = start_pos
        self.end = end_pos
        
    def is_before(self, pos_b):
        """
        Returns true if self comes before Position pos_b without overlapping
        
        Args:
            pos_b = (Position object) another position
        
        Returns: 
            True if Position self is before Position pos_b
            True if they are the same position
            True if chr_b is None
        """
        
        if pos_b == None:
            return True
        
        if self.chr == pos_b.chr:
            return self.end < pos_b.start
            
        return self.chr < pos_b.chr
        
    def overlaps(self, pos_b):
        """
        Returns whether self overlaps Position pos_b
        
        Args: pos_b = (Position object) another position
        
        Returns: 
            (boolean) whether self overlaps with Position pos_b
            False if chr_b is None
        """
        if pos_b == None:
            return False
        
        if self.chr != pos_b.chr:
            return False
        
        start_max = max(self.start, pos_b.start)
        end_min   = min(self.end,   pos_b.end  )
        
        return start_max <= end_min
        
    def starts_before(self, pos_b):
        """
        Returns true if self starts before Position pos_b
        
        Args:
            pos_b = (Position object) another position
        
        Returns: 
            True if Position self starts before Position pos_b
            True if chr_b is None
        """
        
        if pos_b == None:
            return True
            
        if self.chr == pos_b.chr:
            return self.start < pos_b.start
            
        return self.chr < pos_b.chr
    
    def __str__(self):
        return self.chr+":"+str(self.start)+"-"+str(self.end)
        

class Options_list:
    def __init__(self):
    
        #Should lines in the vcf output files be excluded 
        # if they don't have activity?
        # -fv tag sets this to True
        self.filter_vcf = False
        
        #Should lines in the activity (bed) output file be excluded
        # if they don't match a motif?
        # -fa sets this to True
        self.filter_bed = False
        
####-METHODS-####

def get_next_activity(open_activity_file):

    line = open_activity_file.readline().strip()
    
    if line == "":
        return (None, None, None)
    
    line_list = line.strip().split('\t')
    
    #Get position
    chr = line_list[0]
    start_pos = int(line_list[1])
    end_pos = int(line_list[2])
    pos = Position(chr, start_pos, end_pos)
    
    id = line_list[3]
    
    #Get activity scores for samples
    samples_act = []
    for idx in range(4,len(line_list)):
        samples_act.append(float(line_list[idx]))
        
    return (pos, id, samples_act)

def output_activity(open_file, enh_pos, enh_id, matches, options):
    
    #If filtering and there are no matches, print nothing
    if options.filter_bed and len(matches) == 0:
        return
    
    line=enh_pos.chr+'\t'+str(enh_pos.start)+'\t'+str(enh_pos.end)+'\t'+enh_id
    
    for match in matches:
        (name, pos, zscore) = match
        line += '\t'+name+","+str(pos)+","+sf_str(zscore,3)
        
    print(line, file=open_file)

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
    
####-PARSER-####
parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument("-a", "--activity", dest = "activity_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-v", "--vcf", dest = "vcf_dir", required=True)
parser.add_argument("-th", "--threshold", dest = "threshold", 
    required = False, default = 0)
parser.add_argument("-vo", "--vcfout", dest = "vcf_out_dir",
    required=False, default = None)
parser.add_argument("-eo", "--errout", dest = "err_out",
    required=False, default = None)
parser.add_argument("-fa", "--filter_a", action="count", required = False)
parser.add_argument("-fv", "--filter_vcfs", action="count", required = False)

args = parser.parse_args()

# Easier to use argument variables
act_file = args.activity_file
out_file = args.output_file
err_file = args.err_out
vcf_dir = args.vcf_dir
th = float(args.threshold)

options = Options_list()
options.filter_bed = (args.filter_a != None)
options.filter_vcf = (args.filter_vcfs != None)


####-MAIN-####

act_f = open(act_file)
out_f = open(out_file, "w")
err_f = None
if err_file != None:
    err_f = open(err_file, "w")
    print("CHROM\tSTART\tEND\tID\tREFERENCE_ACT\tVARIANTS",file=err_f)

print("CHROM\tSTART\tEND\tID\tVARIANTS",file=out_f)

#Make the vcf output directory if it does not exist
vcf_out_dir_path = None
if args.vcf_out_dir != None and args.vcf_out_dir != vcf_dir:
    vcf_out_dir_path = os.path.join(os.getcwd(),args.vcf_out_dir)
    if not os.path.exists(vcf_out_dir_path):
        os.makedirs(vcf_out_dir_path)
        
headers = act_f.readline().strip().split('\t')
#vcfs will be an array of Vcf objects based on the headers in the activity
# file. The 4 is because the first four columns are not activity scores.
vcfs = [None for idx in range(len(headers)-4)]
#Currently unused
not_found = []

print("Opening vcfs.")
sys.stdout.flush()

#This loop initializes vcfs.
# First four entries are CHROM START END ID
# After that should be samples of the form SAMPLE_OTHERINFO
for idx in range(4,len(headers)):
    #Split on underscore then take the first part
    sample_n = headers[idx].split('_')[0]
    #debug print("SAMPLE: "+sample_n)
    #Find and open the correct vcf based on the sample name
    #getwcd - gets current working directory
    path = os.path.join(os.getcwd(),vcf_dir)
    for vcf_filename in os.listdir(path):
        #Create a Vcf object for the vcf filename that matches the sample name
        vcf_path = os.path.join(path,vcf_filename)
        if os.path.isfile(vcf_path) and sample_n in vcf_filename:
            #Vcf output not currently implemented
            #debug print("Vcf found: <"+vcf_path+">")
            vcf_out_path = None
            if vcf_out_dir_path != None:
                vcf_out_path = os.path.join(vcf_out_dir_path,vcf_filename)
            vcfs[idx-4] = Vcf(sample_n, vcf_path, vcf_out_path)
            break
    if vcfs[idx-4] == None:
        print("No vcf found for sample: "+sample_n+" in directory: "+path)
        not_found.append(sample_n)

sys.stdout.flush()

"""debug
for idx in range(len(vcfs)):
    vcf = vcfs[idx]
    if vcf != None:
        print(str(idx)+" - "+str(vcf.vcf_f))
    else:
        print("None found")
"""

#Enhancer data is of the form:
# (Position object, id, activity list)
#Current enhancer region
(enh_pos, enh_id, enh_act) = get_next_activity(act_f)
print("\tAnalyzing "+enh_pos.chr+"...")
matches = []
# Other enhancers that overlap the current variant
#other_enh = []
# Next enhancer (does not overlap current variant)
#next_enh = None

#debug flag - can set so that only a certain num of matches are calc'd
flag = 0

#Process enhancer regions
while flag < 2 and enh_pos != None:

    #debug print("Enhancer "+enh_id+": "+str(enh_pos))
    #First variant to process (most upstream)
    first_pos = None
    first_idx = -1
    #Overlapping vcfs - vcfs with variants that overlap the enhancer region
    ovl_vcfs = []
    #Scored vcfs for the current enhancer region that have the form:
    # (z-score, information to be output)
    ## Need to decide what info should be output
    scored_vcfs = []
    
    for idx in range(len(vcfs)):
        vcf = vcfs[idx]
        #debug print("\tChecking vcf "+str(vcf))
        #Ensure vcf present and not at end of file
        if vcf == None:
            continue
        #Get next variant from vcf if the variant is upstream of the enhancer
        while vcf.var_pos != None and vcf.var_pos.is_before(enh_pos):
            #debug print("\t\tVariant "+str(vcf.var_pos)+" found to be upstream")
            vcf.next_variant(options)
        
        #debug print("\t\tVar "+str(vcf.var_pos)+" is not upstream")
        
        if vcf.var_pos != None and vcf.var_pos.overlaps(enh_pos):
            #debug print("\t\tMatch found at idx:"+str(idx))
            ovl_vcfs.append(vcf)
            if vcf.var_pos.starts_before(first_pos):
                #debug print(str(vcf)+" "+str(vcf.var_pos)+" ("+str(idx)+
                #    ") comes before "+str(first_pos))
                first_pos = vcf.var_pos
                first_idx = idx
    
    if first_pos != None:
        var_act = enh_act[first_idx]
        #debug print("**First overlapping variant found: "+str(first_pos)+
        #    " (idx = "+str(first_idx)+") with activity: "+str(var_act))
        #Generate lists of 'reference' and samples
        ref_samples_act = []
        #Need to process all overlapping variants together before advancing
        # This will have the form:
        var_samples_idx = []
        for idx in range(len(vcfs)):
            #A 'reference' sample does not have a variant at that position
            # if vcfs[idx] is None then no vcf file could be found so skip
            # if vcfs[idx].var_pos is None then the vcf file has reached its
            #  end and thus it does not have a variant overlapping the
            #  first_pos variant.
            if vcfs[idx] != None and not first_pos.overlaps(vcfs[idx].var_pos):
                #Add the enhancer activity for that sample
                #debug print("\tReference sample includes "+str(vcfs[idx])+
                #    " act = "+str(enh_act[idx]))
                ref_samples_act.append( enh_act[idx] )
            elif vcfs[idx] != None:
                #debug print("\t"+str(vcfs[idx])+" is part of the variant group")
                var_samples_idx.append( idx )
            
        ref_std = 0
        if len(ref_samples_act) > 0:
            ref_std = std(ref_samples_act)
            ref_avg = mean(ref_samples_act)
            #debug print("Std = "+str(ref_std))
            #debug print("Mean = "+str(ref_avg))
        
        if ref_std != 0:
            for idx in var_samples_idx:
                #debug print("Activity = "+str(enh_act[idx]))
                zscore = ( enh_act[idx] - ref_avg ) / ref_std
                #Matches should have name, Position object, and zscore
                if (abs(zscore) > th):
                    matches.append((str(vcfs[idx]) ,vcfs[idx].var_pos ,zscore))
                    vcfs[idx].matches.append((enh_id, enh_pos, zscore))
                #debug print("Z-score for enhancer "+enh_id+" - "+str(enh_pos)+
                #    " and sample "+str(vcfs[idx])+" = "+str(zscore))
        elif err_f != None:
            #Z-score not calculated (insufficient samples or 0 variance)
            err_line = enh_pos.chr+'\t'+str(enh_pos.start)+'\t'
            err_line+= str(enh_pos.end)+'\t'+enh_id+'\t'
            
            err_act = ""
            prev_act = None
            num_times = 0
            for activity in ref_samples_act:
                if prev_act == None:
                    err_act += str(activity)
                    prev_act = activity
                if activity == prev_act:
                    num_times += 1
                else:
                    if num_times > 1:
                        err_act += "x"+str(num_times)
                    err_act += ","+str(activity)
                    
            if prev_act == None:
                err_act = "None"
            else:
                if num_times > 1:
                    err_act += "x"+str(num_times)
            err_line += err_act+"\t"
            
            err_vars = ""
            for idx in var_samples_idx:
                if err_vars != "":
                    err_vars += ","
                err_vars += str(vcfs[idx])
            if err_vars == "":
                err_vars = "None"
            err_line += err_vars
            print(err_line, file=err_f)
            
            
            
        
        for idx in var_samples_idx:
            #debug print("Incrementing "+str(vcfs[idx]))
            vcfs[idx].next_variant(options)

        
    #If no overlapping variant was found, go to next enhancer peak
    else:
        output_activity(out_f, enh_pos, enh_id, matches, options)
        prev_chr = enh_pos.chr
        #debug if len(matches) > 0:
        #    flag += 1
        matches = []
        (enh_pos, enh_id, enh_act) = get_next_activity(act_f)
        if enh_pos != None and prev_chr != enh_pos.chr:
            print("\tAnalyzing "+enh_pos.chr+"...")
    
    sys.stdout.flush()
    
print("Done")
    
act_f.close()
out_f.close()
if err_f != None:
    err_f.close()
for vcf in vcfs:
    if vcf != None:
        vcf.close()