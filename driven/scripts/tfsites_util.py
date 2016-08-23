#!/usr/bin/env python
"""
Last modified: 08/23/2016
Authors: Ethan Pfeifer & Jared Andrews

Calculate motif thresholds for a given motif file and print to a new file.

Usage:
    tfsites_util.py -m <motifs.txt> -o <output.txt>

Args:
    -m (required) = Input filename of a file containing PWMs.

    -o (required) = Output filename.

    -bp (optional) <baselines.txt> = A file containing a single line with tab
        delineated values for baseline probabilities for A, C, G, T (in order).
        Probabilities should all be positive and should sum to 1. If none is
        provided then all are assumed to be equally likely (all are 0.25).

    -pc (optional) <0.1> = Pseudocounts value to be added to all positions of 
        the motif frequency matrix before calculating the probability matrix. 
        
    -th (optional) <0.0> = Default threshold value. This is used if the 
        calculated threshold is lower.
        Ex: default_th = 0.0, biopython calculates threshold needed for
        a given false positive rate is -1.23, threshold printed will be
        0.0.

    -fpr (optional) <0.05> = Acceptable false positive rate for defining 
        thresholds for each motif. 

    -pe (optional) <4> = Precision exponent used by for threshhold calculations. 
        Making this greater than 5 may result in extremely slow run times. Using 
        a lower number will result in faster (but potentially innacurate) 
        calculations.

    -ow (optional flag) = OverWrite: If present, thresholds already present in
        the input file will be replaced in the output file.
"""

from Bio import motifs
import sys
import argparse
parser = argparse.ArgumentParser(usage=__doc__)


def get_baseline_probs(baseline_f):
    """
    Read in baseline probabilities from a file name.
    
    Args:
        baseline_f a file containing a probability array of the form:
        [ PrA PrC PrT PrG ]
        Where PrA + PrC + PrT + PrG = 1 (and all are positive and non-zero)
    
    Returns:
        Array with probabilities as a float array of the form:
        [ PrA, PrC, PrT, PrG ]  
    """
    
    #Default baseline probability numbers (assumes all are equally likely)
    bp_array = [0.25, 0.25, 0.25, 0.25]
    
    with open(baseline_f) as f:
        try:
            for line in f:
                #remove commas, brackets, and whitespace on far left and right
                line = line.strip().replace('[','').replace(']','').replace(',','')
                if line != "":
                    line = line.split()
                    for idx in range(4):
                        bp_array[idx] = float(line[idx])
                    return bp_array
        except ValueError: 
            print ("**ERROR** Baseline probability file incorrectly formatted.\n"+
                "File should contain only [ PrA PrC PrT PrG ] \n"+
                "Where PrA + PrC + PrT + PrG = 1 (and all are positive and non-zero)\n"+
                "Continuing with:")     
            return bp_array
            
        print("**ERROR** Empty file found.\n"+
            "File should contain only [ PrA PrC PrT PrG ] \n"+
            "Where PrA + PrC + PrT + PrG = 1 (and all are positive and non-zero)\n"+
            "Continuing with:")
        return bp_array


def output_motifs(input_f, output, default_th, overwrite, thresholds_list):
    """
    Read in and calculate probability matrices from a frequency matrix file.
    
    Args:
        input = filename of the frequency matrix list
        output = filename to print to
        default_th = default threshold value. This is used if the calculated
            threshold is lower. This value may be None.
            Ex: default_th = 0.0, biopython calculates threshold needed for
            a given false positive rate is -1.23, threshold printed will be
            0.0
        overwrite = should thresholds already written in the file be
            replaced
        thresholds_list = list of thresholds calculated by biopython
        
    
    Returns:
        Motif sequences with names as a list of 
        (id, name, threshold, weight matrix) triples (tuples) from the motif
        (.txt) file. Also returns maximum length motif.
    """

    output_f = open(output,"w")
    idx = 0
             
    # Open provided motif file
    with open(input_f) as f:


        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n 
        # arrays for C, G, T - each with same format as A 
        id = "No id found"
        name = "No name found"
        # Index for line in the motif matrix.
        i = 0
        # Position frequency matrix (read in then reprinted unchanged)
        pfm = ["","","",""]
        given_thresh = None
        
        # Iterate through file line by line.
        for line in f:
            
            #First line contains id and name
            if i == 0:
                line = line.strip().split()
                id = line[0]
                name = line[1]
                #Read in listed threshold
                if (len(line) > 2):
                    given_thresh = float(line[2])
                       
            #Order of position weight matrices is A,C,G,T.
            elif i < 5:       
                pfm[i-1] = line
                
            #Output motif and continue (there are 2 newlines between motifs)
            else:
                if not overwrite and given_thresh != None:
                    th = given_thresh
                elif default_th != None and default_th > thresholds_list[idx]:
                    th = default_th
                else:
                    th = thresholds_list[idx]
                #print the motif back out to the output file
                out_line = id +"\t"+ name +"\t"+ str(th) +"\n"
                out_line += pfm[0] + pfm[1] + pfm[2] + pfm[3]
                print(out_line, file=output_f)
                idx += 1
                i = -1
                given_thresh = None
                pfm = ["","","",""]
            i += 1
    
    if i >= 5:
        if not overwrite and given_thresh != None:
                    th = given_thresh
        elif default_th != None and default_th > thresholds_list[idx]:
            th = default_th
        else:
            th = thresholds_list[idx]
        #print the motif back out to the output file
        out_line = name +"\t"+ id +"\t"+ str(th) +"\n"
        out_line += pfm[0] + pfm[1] + pfm[2] + pfm[3]
        print(out_line, file=output_f)
    
    output_f.close()
    return
 
 
# Create arguments and options
parser.add_argument("-m", "--motif", dest = "motif_file", required=True)
parser.add_argument("-o", "--outfile", dest = "motif_outfile", required=True)
parser.add_argument("-bp", "--baseline", dest = "baseline_file",
    required = False, default = None)
parser.add_argument("-pc", "--pseudocounts", dest = "pseudocounts",
    required = False, default = 0.1)
parser.add_argument("-th", "--threshold", dest = "threshold", 
    required = False, default = None)
parser.add_argument("-fpr", "--falsepos", dest = "false_pos_rate",
    required = False, default = 0.05)
parser.add_argument("-pe", "--precision", dest = "precision_exp",
    required = False, default = 4)
parser.add_argument("-ow", "--overwrite", action="count", required = False)

# Parse arguments and options  
args = parser.parse_args()

if (args.baseline_file == None):
    bp = [0.25, 0.25, 0.25, 0.25]
else:
    bp = get_baseline_probs(args.baseline_file)
pc = args.pseudocounts
if args.threshold != None:
    d_th = float(args.threshold)
else:
    d_th = None
ow = (args.overwrite != None)
fpr = float(args.false_pos_rate)
pe = int(args.precision_exp)
if pe > 5:
    print("Warning: high precision exponent (-pe) may cause drastic slowing "+
    "or memory errors")
if pe <= 0:
    pe = 1
    print("Precision exponent (-pe) too low, set to "+str(pe))
    
thresholds = []
background = {'A':bp[0],'C':bp[1],'T':bp[2],'G':bp[3]}
print("Baseline nucleotide frequencies:\n"+str(background))


print("Calculating thresholds. This could take a while.")
sys.stdout.flush()
idx = 0
exponent = 1

# Calculate thresholds using biopython
fh = open(args.motif_file)
for m in motifs.parse(fh, "jaspar"):
    pwm = m.counts.normalize(pseudocounts = pc)
    pssm = pwm.log_odds(background)
    #Precision argument of 4 was recommended by biopython's documentation
    distribution = pssm.distribution(background=background, precision = 10**pe)
    m_th = distribution.threshold_fpr(fpr)
    thresholds.append(m_th)
    # print progress
    idx += 1
    if (idx >= 10**exponent):
        print( str(idx) + " thresholds calculated..." )
        exponent += 1
        sys.stdout.flush()

print("Total motifs read: "+str(len(thresholds)))
        
print("Outputing thresholds")
# Output calculated motifs
output_motifs(args.motif_file, args.motif_outfile, d_th, ow, thresholds)

print("Done")
