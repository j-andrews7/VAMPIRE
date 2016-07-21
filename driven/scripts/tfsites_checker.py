#!/usr/bin/env python
"""
For a given vcf file (1-based), get reference sequence in which
 it lies and determine if it matches motifs given in a motif file.

Usage: tfsites_checker.py -i <input.vcf> -r <reference.fa>
 -m <motif.txt> -o <output.txt>

Args:
  -i (required) <input.vcf> = Name of variant file to process. Variants must
    be sorted (from low to high) within a chromosome.
  -r (required) <reference.fa> = Name of reference sequence 
    file to get surrounding bases from variant. 3 by default. 
  -m (required) <motif.txt> = Tab-delimited key file containing
    a frequency matrix with each row corresponding to a base
    and each column corresponding to a position
  -o (required) <output.vcf> = Name of output file to be created.

  -pc (optional) <0.1> = Pseudocounts value to be added to all positions of 
    the motif frequency matrix before calculating the probability matrix. 
    Default value is 0.1
  -th (optional) <0> = Motifs are considered a match if they score above a
    given threshold. This is the default threshold (used if no threshold is
    specified by motif file). Default value is 0 (may be a float or int).
  -bp (optional) <baselines.txt> = A file containing a single line with tab
    delineated values for baseline probabilities for A, C, G, T (in order).
    Probabilities should all be positive and should sum to 1. If none is
    provided then all are assumed to be equally likely (all are 0.25)
  -fm (optional flag) = If -fm (filter with motifs) is included, variants
    that do not match any motif will not be included in the output (-o) file.
  -ci (optional) <ChIP.bed> = A bed file containing tab delineated columns
    of the form:
    chr startpos stoppos TF1;TF2;TF3 (etc)
    This must have the start positions in order within each chromosome and 
    must be grouped by chromosome.
  -co (optional) <chip_out.bed> = Name of output bed file to be created.
    A new column will be added with motifs that computationally match each
    peak.
  -fp (optional flag) = If -fk (filter with peaks) is included, ChIP peaks
    that do not match any motif will not be included in the output (-co).
  -sk (optional flag) = If -sk (sorted by karyotype) is included, then the
    program knows that 
"""

import sys
# Used for command line options and such
import argparse
from pyfaidx import Fasta
parser = argparse.ArgumentParser(usage=__doc__)
from math import log2
from math import ceil
from math import log10
from Bio import motifs

####-Classes-####
class Options_list:
    def __init__(self):
        #Should lines in the vcf output file be excluded 
        # if they don't match a motif?
        # -fm tag sets this to True
        self.filter_vcf = False
        #Should lines in the chip (bed) output file be excluded
        # if they don't match a motif?
        # -fp sets this to True
        self.filter_bed = False
        #Is a ChIP peak bed file present?
        # -ci <chip_file.bed> will make this True
        self.chip_present = False
        #Are the input files' chr sorted lexicographically (or by karyotype order)?
        # Input vcf file and chip bed file must be sorted the same way
        # -sk sets this to True
        self.sorted_lex = False

####-Functions-####
def get_surrounding_seq(chromo, var_pos, ref_l, wing_l, fas):
    """ Return sequence containing variant base + specified number
        of bases on each side from reference sequence file.

    Args:
        chromo = Chromosome of variant, e.g. "chr19" or "19" (must be a string)
        var_pos = Integer position of variant start within chromosome.
        ref_l = length of reference sequence. Will be 1 for SNP/insertions but
            greater than 1 for deletions (e.g. deletion of ACTG to G -> ref_l is 4)
        wing_l = Integer number of bases on each side of variant to return (wing
            length) a s full sequence.
        fa_ind = indexed fasta file

    Returns:
        ref_seq = Sequence (string) containing the variant base + specified
        number of bases on each side from reference sequence file.
    """
    
    """#debug 
        print("\tSearching indexed fasta chromosome "+str(chromo)+" at pos "+
            str(var_pos - wing_l - 1)+":"+str(var_pos + wing_l + ref_l - 1))"""
    
    #fas['chr1'][0:1] returns the first base (just one)
    #is either  0 indexed and max (last) base is not returned
    #   or      1 indexed and min (first) base isn't returned
    #A 'sequence' object is returned by fas, but converted by a string for return
    ref_seq = fas[chromo][var_pos - wing_l - 1: var_pos + wing_l + ref_l - 1]

    #debug print("\tSequence: "+str(ref_seq))
    
    return str(ref_seq)
    
def get_reverse_complement(sequence):
    """Returns the reverse complement of a sequence

    Args: sequence = a string containing the original sequence (ACGT only)

    Returns: a string containing the reversed and complementary sequence
    """  
    rev_comp = ""
    
    for idx in range(len(sequence)-1,-1,-1):
        base = 'N'
        
        if sequence[idx].upper() == 'A':
            base = 'T'
        elif sequence[idx].upper() == 'C':
            base = 'G'
        elif sequence[idx].upper() == 'G':
            base = 'C'
        elif sequence[idx].upper() == 'T':
            base = 'A'
        
        rev_comp += base
    
    return rev_comp
        
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
        
def get_baseline_probs(baseline_f):
    """
    Read in baseline probabilities from a file name.
    
    Args:
        baseline_f a file containing a probability array of the form:
        [ PrA PrC PrT PrG ]
        Where PrA + PrC + PrT + PrG = 1 (and all are positive and non-zero)
    
    Returns:
        Array with probabilities as a float array of the form:
        [ PrA, PrC, PrT, PrG ]  """
    
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
             
        
def get_motifs(motif_f, pc, default_th, base_pr):
    """
    Read in and calculate probability matrices from a frequency matrix file.
    
    Args:
        motif_f = filename of the frequency matrix
        pc = pseudocount value (small, greater than zero) to be added to all
            positions in matrix
        default_th = default threshold value used if none is listed in the
            motif file (should be third part of header if present). If None is given,
            default threshold will be calculated for motifs.
        base_pr = probabilities as a float array of the form:
            [ PrA, PrC, PrT, PrG ]
            ** Currently unused by this method but could be used to return a pssm
            instead of a pwm
    
    Returns:
        Motif sequences with names as a list of 
        (id, name, threshold, weight matrix) triples (tuples) from the motif
        (.txt) file. Also returns maximum length motif."""

    motif_list = []
    max_length = 0
             
    # Open provided motif file
    with open(motif_f) as f:


        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n 
        # arrays for C, G, T - each with same format as A 
        id = "No id found"
        name = "No name found"
        #index for which line you are looking at in the motif matrix
        i = 0
        matrix = [[] for x in range(4)]
        
        # Iterate through file line by line.
        for line in f:
        
            #remove brackets and split on whitespace
            line = line.replace('[','').replace(']','').strip().split()
            
            #First line contains id and name
            if i == 0:
                id = line[0]
                name = line[1]
                #If a threshold is given, use it for this matrix. Else, use default.
                if (len(line) > 2):
                    thresh = float(line[2])
                else:
                    thresh = default_th
                
            #Order of position weight matrices is A,C,G,T. Remove brackets etc.
            elif i < 5:       
                matrix[i-1] = line
                if len(matrix[i-1]) > max_length:
                    max_length = len(matrix[i-1])
                
            #add motif to list and continue (there are 2 newlines between motifs)
            else:
                #calculate weights for each position such that each column sums to 1
                p_matrix = calculate_probabilities(matrix, pc)
                tup = (id, name, thresh, p_matrix)
                motif_list.append(tup)
                i = -1
                matrix = [[] for x in range(4)]
            i += 1
            
    
    return (motif_list, max_length)

def calculate_probabilities (matrix, pc):
    """
    Return motif probabilities
    Each column of the input matrix has the number of times that base appears
    at that position in that motif. The occurrences are divided by the total
    number of occurrences to give the frequency of each base.
    Pseudocounts are added in to avoid having a zero value in matrices with
    small sample sizes
    """
    #for clarity
    a = matrix[0]
    c = matrix[1]
    t = matrix[2]
    g = matrix[3]
    
    #initialize output matrix of same size as input matrix
    new_m = [[0 for y in range(len(a))] for x in range(4)]
    
    #divide each position by the sum of its column
    for pos in range(len(a)):
        total = float(a[pos]) + float(c[pos]) + float(t[pos]) + float(g[pos])
        total += 4 * pc
        for base in range(4):
            new_m[base][pos] = ( float(matrix[base][pos]) + pc ) / total
    
    return new_m

    
def match_motifs(motifs,b_prob,ref_seq,var_seq,wing_l):
    """
    Takes a reference and variant sequence, then checks if they match the motif
    list. Outputs match score for both ref seq and var seq in any case where
    either matches above the given threshold.
    
    Args:
        motifs = Motif sequences with names as a list of the form:
            (id, name, threshold, weight matrix) tuples from the motif 
            (.txt) file. Also returns maximum length motif.
        b_prob = sequence of baseline probabilities of each base, in order (A, C,
            G, T). Probabilities should sum to 1.
        ref_seq = Reference sequence represented as a string (ACGT only)
        var_seq = Variant sequence represented as a string (ACGT only)
        wing_l = Integer length of sequence of bases flanking the variant
        
    Returns: A list of tuples of the form:
        (motif id, motif variant score, motif reference score, ChIP match)
            Where: motif id is a string 
            motif score is a float
            ChIP match is a boolean (False here until matched later)
    """
    
    #motif match array with the form (id, name, p_matrix, max_score, match_seq)
    scored   = score_motifs(motifs,b_prob,var_seq,wing_l)
    r_scored = score_motifs(motifs,b_prob,ref_seq,wing_l)
    
    matches = []
    
    #generate list of motifs that matched var seq to compare ref seq matches
    for idx in range(len(scored)):
        ( id,  name,  th,  max_score,  match_seq)  = scored[idx]
        (rid, rname, rth, rmax_score, rmatch_seq) = r_scored[idx]
        if (max_score >= th or rmax_score >= rth):
            #debug
            if (id != rid or name != rname or th != rth):
                print("***ERROR*** matching motifs to varseq and refseq desynced\n"+
                    id +" != "+rid+" or "+name+" != "+rname+" or "+th+" != "+rth)
            #tup = (id, name, max_score, match_seq, rmax_score, rmatch_seq)
            tup = (name, max_score, rmax_score, False)
            matches.append(tup)
    
    #output variant seq matches vs ref seq matches
    """for m in matches:
        line = "Variant matched motif "+m[1]+" ("+m[0]+") with score "
        line+= str(m[2])[:6]+" at seq: "+m[3]+"\nReference matched motif "+m[1]
        line+= "           with score "+str(m[4])[:6]+" at seq: "+m[5]
        #debug line+= "\n\tRefseq: "+ref_seq+"\n\tVarseq: "+var_seq
        print(line,file=output_f)"""
        
    return matches
    
def score_motifs(motifs, baseline_p, sequence, wing_l):
    """ 
    Calculate if any motifs in the motif list match the given sequence. Requires
    that no motif have a length of 0. 
    
    Args:
        motifs = Motif sequences with names as a list of (id, name, weight matrix)
            triples (tuples) from the motif (.txt) file. Also returns maximum length
            motif.
        baseline_p = sequence of baseline probabilities of each base, in order (A, C,
            G, T). Probabilities should sum to 1.
        sequence = array of characters (strings 'A','C','G', or 'T')
        wing_l = length of sequence of bases flanking the variant
    
    Returns:
        matches = array of matches between motifs and sequences of the form
            (id, name, thresh, max_score, match_seq)
    
    Concern: Longer motifs matches will generate higher scores
    """
    #list of matches between motifs and sequence
    matches = []
    #print("Running match motifs!")
    
    for tuple in motifs:
        (id, name, thresh, p_matrix) = tuple
        #highest match score for this motif (based on starting position)
        max_score = float("-inf")
        #sequence that matched motif best
        match_seq = ""    
        
        #trim flanking bases to one less than motif length
        motif_l = len(p_matrix[0])
        #number of bases to trim
        trim_amount = wing_l-motif_l+1
        #sequence with flanking bases short enough that motif overlaps variant
        trim_seq = sequence
        
        if trim_amount != 0:
            trim_seq = sequence[trim_amount:-trim_amount]
            
        """#debug trimming
        if trim_amount != 0:  
            print("Sequence trimmed from:\'"+sequence+"\' to \'"+trim_seq+"\' for "+
                str(motif_l)+" length motif")
        else:
            print("Sequence '"+sequence+"\' not trimmed for "+str(motif_l)+
                " length motif")"""
        
        #iterate and check match to all positions where motif overlaps with variant
        for pos in range(motif_l):
            #check match starting at position (score_motif will stop after the length
            # of the motif)
            pos_score = score_motif(p_matrix, baseline_p, trim_seq[pos:])
            if pos_score > max_score:
                max_score = pos_score
                match_seq = trim_seq[pos:pos+motif_l]
        
        #debug print("Max match score:"+str(max_score)+" for motif "+name+" and sequence "+match_seq+".")
        tup = (id, name, thresh, max_score, match_seq)
        matches.append(tup)
        
        
    return matches
        
def score_motif(p_matrix, baseline_p, sequence):
    """
    Calculate match between given probability matrix and given sequence starting
        at the beginning of the given sequence and ending after the length of the
        motif.
    
    Args:
        p_matrix = probability matrix where rows are bases (A,C,G,T) and columns
            are probabilities (that sum to 1) at a given position. All rows should
            be the same length.
        baseline_p = sequence of baseline probabilities of each base, in order (A, C,
            G, T). Probabilities should sum to 1.
        sequence = array of characters (strings 'A','C','T', or 'G')
    
    Returns:
        Match score between probability matrix and sequence
    
    Concern: Longer motifs matches will generate higher scores
    """
    
    score = 0
    
    for pos in range(len(p_matrix[0])):
        #Match base
        if   (sequence[pos] == 'A' or sequence[pos] == 'a'):  
            # natural log of likelihood ratio
            
            score += score_base(p_matrix[0][pos], baseline_p[0])
            
        elif (sequence[pos] == 'C' or sequence[pos] == 'c'):
            score += score_base(p_matrix[1][pos], baseline_p[1])
        elif (sequence[pos] == 'G' or sequence[pos] == 'g'):
            score += score_base(p_matrix[2][pos], baseline_p[2])
        elif (sequence[pos] == 'T' or sequence[pos] == 't'):
            score += score_base(p_matrix[3][pos], baseline_p[3]) 
        else:
            #print("non ACGT base\'"+sequence[pos]+
            #  "\' found in sequence \'"+sequence+"\'.")
            sys.stdout.flush()
     
    """#debug calculations for a score match
        if score > 0:
            print("*Match score of "+str(score)+" for "+sequence)
            calcs = "*Calculations:"
            score = 0
            for pos in range(len(p_matrix[0])):
                if pos != 0:
                    calcs += " + "
                if pos%5 == 0:
                    calcs += "\n"
                if   (sequence[pos] == 'A' or sequence[pos] == 'a'):  
                    score += score_base(p_matrix[0][pos], baseline_p[0])
                    calcs += str(score_base(p_matrix[0][pos], baseline_p[0]))[:5]
                elif (sequence[pos] == 'C' or sequence[pos] == 'c'):
                    score += score_base(p_matrix[1][pos], baseline_p[1])
                    calcs += str(score_base(p_matrix[1][pos], baseline_p[1]))[:5]
                elif (sequence[pos] == 'G' or sequence[pos] == 'g'):
                    score += score_base(p_matrix[2][pos], baseline_p[2])
                    calcs += str(score_base(p_matrix[2][pos], baseline_p[2]))[:5]
                elif (sequence[pos] == 'T' or sequence[pos] == 't'):
                    score += score_base(p_matrix[3][pos], baseline_p[3])
                    calcs += str(score_base(p_matrix[3][pos], baseline_p[3]))[:5]
            calcs += "\n = "+str(score)[:5]
            print(calcs)
        """
        
    return score
 
def score_base(probability, baseline_p):
    if (probability == 0):
        #should mathematically be negative infinity
        #this should not occur with pseudocounts added in
        return -100
    else:
        return log2(probability / baseline_p)

        
def get_next_var(opened_file):
    """
    Reads in the next line of the vcf and returns the next variant's information
    
    Args: opened_file = an already open input .vcf file
    
    Returns: a tuple with the following information (in order) or None
        chromosome number as a string e.g. "chr1",
        position of the variant as an integer
        reference sequence
        variant sequence
        """
    
    line = opened_file.readline()
    
    #skip info lines-
    while line.startswith("##"):
        line = opened_file.readline()
    
    line = line.strip()
    
    #input file is empty
    if line == "":
        return None
    
    line_list = line.split("\t")
    
    chr = line_list[0]
    pos = int(line_list[1])
    ref = line_list[3]
    var = line_list[4]

    return (chr, pos, ref, var, line)
            
def update_vcf(line, matches, output_f, options):
    """
    Updates the output file with the output in the correct variant call format
    
    Args:
        line = original vcf line (input file line)
        matches = list of motif matches as tuples of the form:
            (id, variant score, reference score, ChIP match)
        output_f = output vcf
        options = options list
    
    Returns: Nothing (updates output_f instead of returning)
    """
    
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
    
    
    for match in matches:
        if chips != "":
            names += ","
            varscores += ","
            refscores += ","
            chips += ","
        names += match[0]
        varscores += sf_str(match[1],3)
        refscores += sf_str(match[2],3)
        if match[3]:
            chips += "Y"
        else:
            chips += "N"
    
    #If there are no matches, print the line unchanged or filter it out (return
    #without printing)
    if len(matches) == 0:
        if options.filter_vcf == False:
            print(line, file=output_f)
        return
    
    outline = ""
    idx = 0
    
    for col in columns:
        if outline != "":
            outline += "\t"
        if idx == 7:
            outline += "MOTIFN="+names+";MOTIFV="+varscores+";MOTIFR="+refscores
            if (options.chip_present):
                outline += ";MOTIFC="+chips
            if col != '.':
                outline += ";" + col
        else:
            outline += col
        idx += 1
    
    if idx < 7:
        print("**Error** VCF formatted incorrectly. Less than 8 columns found:\n"
            +line)
        #Add output at the end anyway
        outline += "\tMOTIFN="+names+";MOTIFV="+varscores+";MOTIFR="+refscores
        if (options.chip_present):
            outline += ";MOTIFC="+chips
    
    print (outline, file=output_f)
    
    return

    
def match_peaks(chr, pos, peaks, chip_f, matches, output_f, options):
    """
    Returns an array of peaks that match the current chromosome and position.
    Updates the coutput_f if one is present.
    
    Args:
        p_chr = (string) previous chromosome. Needed to know if a new chromosome
            is being entered.
        chr = (string) chromosome. Chromosomes 1-22, X, and Y are expected.
        pos = current position of the variant.
        peaks = buffer of peaks. They should all be upstream of or overlapping the
            variant at chr and pos. Peaks is an array of tuples of the form:
            (chr, start pos, end pos, array of ChIP tfs, motif match array)
            tfs = transcription factors
        chip_f = input ChIP bed file to be read from. This must have the start
            positions in order within each chromosome and must be grouped by
            chromosome.
        matches = list of motif matches as tuples of the form:
            (name, variant score, reference score, ChIP match)
        output_f = ChIP output bed file to be printed to.
        options = Options_list object
    
    Returns: Peak buffer tuple of the form ( overlapping peak array, next peak )
        Array of peaks that overlap the current chromosome and position
        Next peak (because you have to over-read to make sure you don't miss any)
        match_peaks also updates coutput_f.
    """
    
    if chip_f == None:
        return ([],matches)
    
    #Get rid of peaks that are upstream of the current chromosome
    idx = 0
    while idx < len(peaks) and chr_less(peaks[idx][0], chr, options):
        #If the chromosome doesn't match, output the line and keep searching
        print_peak(peaks[idx], output_f, options)
        idx += 1

    #peak at idx will be included and the rest will be removed
    peaks = peaks[idx:]
    
    #If previous peaks were not from correct chromosome, get there
    if (len(peaks) == 0):
        new_peak = get_peak_at(chr, pos, chip_f, output_f, options)
        #If end of file is reached
        if (new_peak == None):
            return ([],matches)
        else:
            peaks.append(new_peak)
    
    idx = 0
    
    #Read through bed file 
    while True:

        #If more peaks are needed
        if idx == len(peaks):
            n_peak = get_next_peak( chip_f ) 
            
            #If end of bed file is reached, then just return current list
            if n_peak == None:
                return (peaks,matches)
                
            peaks.append(n_peak)

        #Current peak (chromosome, start pos, end pos, transcription factor array,
        # matrix match array)
        (pchr, psta, pend, ptfs, pmms) = peaks[idx]
        
        #If next chromosome is reached in bed file
        if pchr != chr:
            break
        
        if psta <= pos:
            if pend >= pos:
                motif_idx = 0
                for motif_idx in range(len(matches)):
                    (m_name, mv, mr, m_chip) = matches[motif_idx]
                    pmms.append(matches[motif_idx])
                    for trans_factor in ptfs:
                        #If the transcription factor (chip peak) name is the same as
                        # the matched motif name, note that there is a chip match
                        if trans_factor == m_name:
                            #Motif match is verified by ChIP data
                            matches[motif_idx] = (m_name, mv, mr, True)
                #Save with new value for pmms
                peaks[idx] = (pchr, psta, pend, ptfs, pmms)
            #Otherwise both are before pos, so remove that peak and continue
            #This should only ever happen when idx is 0... but still
            else:
                print_peak(peaks[idx], output_f, options)
                peaks = peaks[0:idx]+peaks[idx+1:]
                idx -= 1
        #Otherwise peak start is after the variant position, so stop
        else:
            break;
        idx += 1
    
    return (peaks,matches)
         
def get_next_peak( opened_file ):
    """
    Reads in the next line of the bed and returns the next peak's information
    
    Args: opened_file = an already open input .bed file
    
    Returns: a tuple with the following information (in order) or None
        chromosome number as a string e.g. "chr1",
        position of the start of the peak
        position of the end of the peak
        array list containing transcription factors which bind in that area
        array list containing motif matches (empty)
    """
    
    #print("Entering get_next_peak")
    #sys.stdout.flush()
    
    line = opened_file.readline().strip()
    
    #print("Got: <"+line+">")
    sys.stdout.flush()
    
    #input file is empty
    if line == "":
        return None
    
    line_list = line.split('\t')

    chr = line_list[0]
    start = int(line_list[1])
    end = int(line_list[2])
    tf_array = line_list[3].split(';')
    
    return (chr, start, end, tf_array, [])
    
def get_peak_at( chr, pos, chip_f, out_f, options ):
    """
    Get the first peak where the end point of the peak is past the input pos.
    Requires that chip_f is sorted the same way is the vcf input file.
        This file will print all intermediate peaks to the out_f if one is given.
        
    Args:
        chr = the new chromosome to get to
        pos = the new position to get to
        chip_f = an already open input .bed file
        out_f = an already open output file (to be printed to). May be None.
        options = object from Options_list class
        
    Returns:
        The first peak from chip_f that is the same chromosome as chr
        where the end point of the peak is past the input position (pos).
        If the end of file is reached, None is returned.
    """
    
    #print("*Entering get_peak_at")
    #sys.stdout.flush()
    
    #Skip ahead until the correct chromosome
    peak = get_next_peak( chip_f )
    
    while peak != None:
        (p_chr, p_sta, p_end, p_tfa, p_mm) = peak
        #If the chromosome is correct, skip ahead until the correct position
        if p_chr == chr:
            #We have passed the position at the chromosome
            if p_end >= pos:
                #print("get_peak_at returns: "+p_chr+":"+str(p_sta)+"-"+str(p_end))
                #sys.stdout.flush()
                return peak
            else:
                print_peak(peak, out_f, options)
                
        #If the chromosome is too low and there is an outfile, print
        elif chr_less(p_chr, chr, options):
            print_peak(peak, out_f, options)
        #If we have passed the chromosome
        else:
            return peak
            
        peak = get_next_peak( chip_f )
    #print("*get_peak_at returns None")
    #sys.stdout.flush()
    return None
                        
def print_peak( peak, file, options ):
    """
    Prints the peak to the given file (or exits if no file is given)
    
    Args:
        peak = ChIP peak of the form:
            (chr, start, stop, chip tf array, motif match tf array)
            chip tf array is an array of tf names
            motif match tf array is an array of (motif name, vscore, rscore, strand) 
        file = the file to print to. This file should already be opened. If this
            file is None, nothing will happen.
        options = Options_list object
        
    Returns: Nothing
    """
    if file == None:
        return
    
    (chr, start, end, c_array, mm_array) = peak
    
    #If there are no motif matches and filtering is on, do not print this peak
    if options.filter_bed and len(mm_array) == 0:
        return
    
    line = chr +'\t'+ str(start) +'\t'+ str(end)+'\t'
    
    #Generate string of chip transcription factors from the peak
    chip_string = ""
    for tf_name in c_array:
        if chip_string != "":
            chip_string += ";"
        chip_string += tf_name
    
    #Generate string of motif matches that overlap the peak
    motif_string = ""
    for match in mm_array:
        if motif_string != "":
            motif_string += ";"
        (name, vscore, rscore, chip) = match
        motif_string += name+","+sf_str(vscore,3)+","+sf_str(rscore,3)
    
    line += chip_string + "\t" + motif_string
    print(line, file=file)
    
    return
        
def chr_less( chr_left, chr_right, options):
    """
    Returns true if the left chromosome comes before the right or is the same.
    
    Args:
        chr_left = (string) the left chromsome being compared
        chr_right = (string) the right chromosome being compared
        options = Options_list object
    
    Returns: Whether or not the left chromosome comes before the right (boolean).
    """
    
    #True if the file is sorted lexicographically
    # i.e. chr1 < chr11 < chr2 < chrX < chrY
    if options.sorted_lex:
        return (chr_left < chr_right)
        
    #False if the file is sorted numerically
    # i.e. chr1 < chr2 < chr11 < chrX < chrY
    else:
        left = chr_left[3:]
        right = chr_right[3:]
        try:
            l_num = int(left)
            try:
                r_num = int(right)
                return l_num < r_num
            #Right chromosome is a string (chrX, chrY, chrTest, etc)
            except:
                #Left is a number and right is a string
                #Numbers are sorted before strings (chr1 < chrX)
                return True
        #Left number is a string
        except ValueError:
            try:
                r_num = int(left)
                #Left is a string and right is a number
                #Numbers are sorted before strings (chrX !< chr1)
                return False
            #Both are strings, sort lexicographically
            except ValueError:
                return chr_left < chr_right
                
                
    
####-PARSER-####

""" Requires that vcf file have variants sorted by position within chromosomes
"""

# Create arguments and options
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-r", "--ref", dest = "ref_file", required=True)
parser.add_argument("-m", "--motif", dest = "motif_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-ci", "--chip", dest = "chip_file",
    required = False, default = None)
parser.add_argument("-co", "--chipout", dest = "chip_out_file",
    required = False, default = None)
parser.add_argument("-bp", "--baseline", dest = "baseline_file",
    required = False, default = None)
parser.add_argument("-th", "--threshold", dest = "threshold", 
    required = False, default = None)
parser.add_argument("-pc", "--pseudocounts", dest = "pseudocounts",
    required = False, default = 0.1)
parser.add_argument("-fm", "--filter_o", action="count", required = False)
parser.add_argument("-fp", "--filter_co", action="count", required = False)
parser.add_argument("-sk", "--kary_sort", action="count", required = False)

args = parser.parse_args()

# Easier to use argument variables
inp_file = args.input_file
ref_file = args.ref_file
motif_file = args.motif_file
out_file = args.output_file
bp_file = args.baseline_file
chip_file = args.chip_file
chip_out_file = args.chip_out_file
pc = float(args.pseudocounts)

#Options list. Easier to pass in methods or use in code updates.
options = Options_list()
options.filter_vcf = (args.filter_o != None)
options.filter_bed = (args.filter_co != None)
options.sorted_lex = (args.kary_sort == None)

#Output so user can double check options
print(("Input file: {}\nReference file: {}\nMotif file: {}\n"+
    "Output file: {}\nOptional arguments:\n    Pseudocounts value = {}").
    format(
        inp_file, ref_file, motif_file,
        out_file, pc
        ))
    
#Optional artuments to print
opt_args = ""

if args.threshold != None:
    th = float(args.threshold)
    opt_args += "    Default match threshold = "+str(th)+"\n"
else:
    th = 0.0
    opt_args += "    Default match threshold = "+str(th)+"\n"

chip_f = None
coutput_f = None

if (chip_file != None):
    options.chip_present = True
    opt_args += "    ChIP file: "+chip_file+"\n"
    chip_f = open(chip_file)
    if (chip_out_file != None):
        opt_args += "    ChIP output file: "+chip_out_file+"\n"
        coutput_f = open(chip_out_file,"w")
elif (chip_out_file != None):
    opt_args += "No ChIP file given, so no ChIP output file will be created\n"

if (bp_file != None):
    opt_args += "    Baseline probabilities file: "+bp_file+"\n"
if (options.filter_vcf):
    opt_args += "    Filter output vcf for motif matches? Yes\n"
if (options.filter_bed):
    opt_args += "    Filter output ChIP bed for motif matches? Yes\n"
if (not options.sorted_lex):
    opt_args += "    Input vcf and Input ChIP bed are sorted by karyotype\n"

print(opt_args)
    
####-Main-####

if (bp_file == None):
    bp = [0.25, 0.25, 0.25, 0.25]
else:
    print("Reading in baseline probabilities:")
    bp = get_baseline_probs(bp_file)
    print(bp)
sys.stdout.flush()

# Grab motif list from motif file.
print("Creating motif list from " + motif_file)
sys.stdout.flush()
(motifs, max_motif_l) = get_motifs(motif_file, pc, th, bp)
print("Maximum motif length is "+str(max_motif_l)+".")


#Wing length should be one less than the length of the maximum motif so
#motif is only matched against places with overlap with the variant
wing_l = max_motif_l - 1

# Open output file.
output_f = open(out_file,"w")


"""#debug that motifs were calculated correctly
    for motif in motifs:
        print(motif[0] + " " + motif[1]+" "+str(motif[2]))
        sys.stdout.flush()
        for i in range(4):
            line = "[ "
            for p in range ( len( motif[3][i] ) ):
                line+= sf_str(motif[3][i][p],3) + "  "
            print(line+"]")
        print("")
    print("")"""

# Create index file from input fasta for quick searching
print("Creating index from reference sequence for efficient searching...")
print("This will be slow the first time.")
sys.stdout.flush()
fa_ind = Fasta(ref_file)

print("Analyzing variants. This may take a while.\n")
sys.stdout.flush()

# Open VCF file.
with open(inp_file) as vcf:

    line = vcf.readline()
    
    info_needed = True
    info = "##INFO=<ID=MOTIFN,Number=.,Type=String,Description="
    info += "\"Matched motif names\">\n"
    info += "##INFO=<ID=MOTIFV,Number=.,Type=Float,Description="
    info += "\"Variant match scores\">\n"
    info += "##INFO=<ID=MOTIFR,Number=.,Type=Float,Description="
    info += "\"Reference match scores\">\n"
    if (chip_f != None):
        info += "##INFO=<ID=MOTIFC,Number=.,Type=Character,Description="
        info += "\"Motif validated by ChIP (Y/N)\">"
    
    #Skip info lines
    while line.startswith("##"):
    #Print new info lines at the top of the ##INFO section
        if info_needed and line.startswith("##INFO"):
            print(info, file=output_f)
            info_needed = False
        print(line, file=output_f, end="")
        line = vcf.readline()
    
    # Create appropriate header.
    header = line.strip()
    print(header, file=output_f)

    #Previous computed sequence
    p_chr = ""
    p_end = 0
    p_seq = ""
    
    curr_var = get_next_var(vcf)
    
    #Queue of upcoming variants that may be close enough to require modification
    # of the surrounding sequence downstream of the current variant
    #Will contain variants stored as a tuple of the form:
    # (chromosome, start position, reference seq, variant seq, line text)
    next_vars = []
    
    #Queue of ChIP peaks that overlap the current variant
    #Will contain peaks as a tuple of the form
    # (chr, start, stop, chip tf array, motif match tf array)
    # chip tf array is an array of tf names
    # motif match tf array is an array of (motif name, vscore, rscore, strand)
    peak_buffer = []
    
    #Process each variant
    while True:
        
        (chr, pos, ref_bases, var_bases, line) = curr_var
        
        #Get reference sequence surrounding the variant from the reference file
        ref_seq = get_surrounding_seq(chr, pos, len(ref_bases), wing_l, fa_ind)
        surr_seq = ref_seq
        
        #check that reference sequence matches the vcf's sequence for that position
        returned_ref_b = surr_seq[wing_l:-wing_l]
        if ref_bases.upper() != returned_ref_b.upper():
            print("**ERROR**\nVCF reference sequence for "+chr+" pos "+str(pos)+
                ":\n\""+ref_bases+"\" does not match reference file sequence:\n\""+
                returned_ref_b+"\"")
        
        #debug sequence generation
        #print("analyzing "+chr+":"+str(pos)+" "+ref_bases+"->"+var_bases)
        #print(surr_seq)
        
        #check if the last sequence overlapped the current sequence
        if chr == p_chr and p_end >= pos - wing_l :
            #Don't forget to count the overlapping end base
            overlap = p_end - (pos - wing_l) + 1
            """#debug overlap
                print("Calculating overlap from prev seq ending at "+str(p_end)+": \'"+
                    p_seq+"\' and \'"+surr_seq+"\'\n\tat pos "+str(pos)+" with overlap of "
                    +str(overlap))"""
            overlap_seq = p_seq[-overlap:]
            surr_seq = overlap_seq + surr_seq[overlap:]
            #debug overlap print("New sequence is \'"+surr_seq+"\' with overlap of "+str(overlap))
        
        #position of the last base in the surrounding sequence
        end_pos = pos + wing_l + len(ref_bases)
        #debug print("\tEnd pos calculated to be:"+str(end_pos))
        
        #check upstream for variants that would affect surrounding sequence
        #will need to check a maximum of wing_l number of variants because
        # each should leave behind at least one base
        for idx in range(0,wing_l):
            #if the next variant is not in the list, retrieve it from the file
            if idx == len(next_vars):
                next_var = get_next_var(vcf)
                if next_var == None:
                    break
                else:
                    next_vars.append(next_var)    
             
            (n_chr, n_pos, n_ref, n_var, nl) = next_vars[idx]
            
            #check if the next variant overlaps the surrounding sequence of the
            # current variant
            if (n_chr == chr and n_pos <= end_pos):
                """#debug next variants 
                    print("Downstream variant found "+n_chr+":"+str(n_pos)+" "+n_ref+"->"+
                        n_var)"""
                    
                #The offset is the start point of the next variant in the string
                #And it needs to be counted from the end of the wing (is negative)
                offset = n_pos - end_pos
                #debug print("\tOffset is: "+str(offset)+", End pos is: "+str(end_pos))
                #End point of the same variant in the string. Should also be neg.
                if (offset + len(n_ref) < 0):
                    surr_seq = surr_seq[:offset] + n_var + surr_seq[offset+len(n_ref):]
                else:
                    #variant sequence of deletion extends past the end of surr seq
                    surr_seq = surr_seq[:offset] + n_var
                    #end point must be updated (is the endpoint of the deletion)
                    #debug print("\tEnd position was: "+str(end_pos))
                    end_pos += offset + len(n_ref)
                    #debug print("\tNew end position: "+str(end_pos))
                    
                #debug next variants print("Downstream variant added, surr_seq = \'"+surr_seq+"\'")
                
                #in the case of a deletion, more reference sequence is now needed
                bases_needed = (2*wing_l + len(ref_bases)) - len(surr_seq)
                
                #add extra bases and update end position of sequence
                if bases_needed > 0:
                    surr_seq += str(fa_ind[chr][end_pos - 1:end_pos + bases_needed - 1])
                    end_pos += bases_needed
                    """#debug bases_needed
                        print(str(bases_needed)+" reference bases added,"+
                            " surr_seq = \'"+surr_seq+"\'")"""
                    
            #stop checking if next variant isn't in wing sequence of current variant
            else:
                break
                        
        #update to replace current reference bases with reference bases
        surr_seq = surr_seq[:wing_l] + var_bases + surr_seq[wing_l+len(ref_bases):]
        
        #Variant sequence with correct wing sizes 
        #(right wing may have been by a downstream insertion)
        var_seq = surr_seq[:2*wing_l+len(var_bases)]
        
        #debug print("Calculated variant sequence: "+var_seq+", matching motifs...")
        
        #Calculate motif matches to variant sequence
        plusmatch = match_motifs(motifs,bp,ref_seq,var_seq,wing_l)
        
        #Also output matches to reverse complement
        r_refseq = get_reverse_complement(ref_seq)
        r_varseq = get_reverse_complement(var_seq)
        
        #Calculate motif matches to reverse complement
        minusmatch = match_motifs(motifs,bp,r_refseq,r_varseq,wing_l)
        
        matches = plusmatch + minusmatch
        
        #print("calling match_peaks for variant "+chr+":"+str(pos))
        #sys.stdout.flush()
        
        #Update ChIP buffer for current position
        #Update matches array with peak overlap data
        (peak_buffer, matches) = match_peaks(chr, pos, peak_buffer, chip_f,
            matches, coutput_f, options)
        
        """print("match_peaks returned "+str(len(peak_buffer))+" peak(s):")
        for peak in peak_buffer:
            (pchr, psta, pend, ptfs, pmms) = peak
            print(pchr+":"+str(psta)+"-"+str(pend)+" tfs:"+str(len(ptfs))+
                " mms:"+str(len(pmms)))
        print()"""
        
        #Create the correct line in VCF format and print to output_f
        update_vcf(line, matches, output_f, options)
        
        sys.stdout.flush()
        
        #Update previous, next, and current variables
        if (p_chr != chr):
            print("\tAnalyzing "+chr+"...")
        
        p_chr = chr
        #position of the last base included
        p_end = pos + len(ref_bases) - 1
        #p_seq only holds sequence related to bases upstream of the next variant
        p_seq = surr_seq[:wing_l+len(var_bases)]
        
        #pop(0) takes from the end of the queue (front of the sequence)
        if (len(next_vars) > 0):
            curr_var = next_vars.pop(0)
        else:
            break

    #Print remaining peaks
    for peak in peak_buffer:
        print_peak(peak, coutput_f, options)
        

        
"""#debug motif matching  
    #wing = "AAAAAAAAAA"
    #match_motifs(motifs,[0.25,0.25,0.25,0.25],wing+"GTCTGTGGTTT"+wing,len(wing))
    #match_motifs(motifs,[0.25,0.25,0.25,0.25],wing+"CACGTG"+wing,len(wing))"""
 
print("Done")

# Close output files.
output_f.close()
if coutput_f != None:
    coutput_f.close()
    
