#!/usr/bin/env python3
"""
For a given vcf file, compare each variant and its reference sequence to
a set of motifs to determine if any are significantly altered.

Usage:
    python motifs.py -i <input.vcf> -r <reference.fa> -m <motif.txt>
        -o <output.vcf> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process.
    -r (required) <reference.fa>: Name of reference sequence
        file to get surrounding bases from variant.
    -m (required) <motif.txt>: Tab-delimited key file containing a frequency
        matrix with each row corresponding to a base and each column
        corresponding to a position (JASPAR format).
    -o (required) <output.vcf>: Name of output file to be created.

    -pc (optional) <0.1>: Pseudocounts value to be added to all positions of
        the motif frequency matrix before calculating the probability matrix.
    -th (optional) <0>: Motifs are considered a match if they score above a
        given threshold. This is the default threshold (used if no threshold is
        specified by motif file).
    -ws (optional) <50>: Wing size in bp to search for weak homotypic
        matches, co-binding tfs, and GC content.
    -bp (optional) <baselines.txt>: A file containing a single line with tab
        delineated values for baseline probabilities for A, C, G, T (in order).
        Probabilities should all be positive and should sum to 1. If none is
        provided then all are assumed to be equally likely.
    -fm (optional flag): If -fm (filter with motifs) is included, variants
        that do not match any motif will not be included in the output file.
    -ci (optional) <ChIP.bed>: A sorted bed-like file containing tab delineated
        columns of the form:
        chr start end TF1;TF2;TF3...
    -co (optional) <chip_out.bed> = Name of output bed file to be created.
        A new column will be added with motifs that computationally match each
        peak.
    -fp (optional flag): If -fp (filter with peaks) is included, ChIP peaks
        that do not match any motif will not be included in the output (-co).
    -sk (optional flag): Use if sorted by karyotype.
    -fc (optional): filter_chip YYY; can not also be called with fn
    -fn (optional): filter_novel YYY; can not also be called with fc
    -mv (optional) <ws>: merge adjacent variants for same chromosome and sample
        that occur within passed integer distance. default distance = wing size.
"""

import sys
import argparse
import time
import pdb    # necessary for debugger; use pdb.set_trace()
import motif
import sequence
from pyfaidx import Fasta
from Bio import motifs

parser = argparse.ArgumentParser(usage=__doc__)


# TODO - Remove use of this, will only create headaches later.
class Options_list:

    def __init__(self):
        # Should lines in the vcf output file be excluded
        # if they don't match a motif?
        # -fm tag sets this to True
        self.filter_vcf_motif = False
        # Should lines in the vcf output file be excluded
        # if they don't match a ChIP peak?
        # -fc tag sets this to True
        self.filter_vcf_chip = False
        # Should lines in the vcf output file be excluded
        # if they do match a ChIP peak? Allows for printing of only
        # potentially novel variants
        # -fn tag sets this to True
        self.filter_vcf_no = False
        # Is a ChIP peak bed file present?
        # -ci <chip_file.bed> will make this True
        self.chip_present = False


def timeString():
    """ Return time as a string YearMonthDay.Hour.Minute.Second
    """
    return time.strftime('%Y%m%d.%H:%M:%S')


def update_vcf(line, matches, output_fileHandle, options):
    """
    Updates the output file with the output in the correct variant call format

    Args:
        line = original vcf line (input file line)
        matches = list of MotifMatch objects
        output_fileHandle = output vcf file handle
        options = options list

    Returns: Nothing (updates output_fileHandle instead of returning)
    """

    # First 8 columns should always be:
    # CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    line = line.strip()
    columns = line.split('\t')

    # ID=MOTIFN,Type=String,Description="Matched motif names"
    # ID=MOTIFV,Type=Float,Description="Motif variant match scores"
    # ID=MOTIFR,Type=Float,Description="Motif reference match scores"
    # ID=MOTIFC,Type=Character,Description="Motif validated by ChIP (Y/N)"
    names = ""
    varscores = ""
    refscores = ""
    chips = ""
    varht = ""
    refht = ""
    vargc = ""
    refgc = ""

    to_be_printed = 0

    for match in matches:
        # First, check if the match needs to be filtered
        if options.filter_vcf_chip and not match.chip_match:
            continue
        elif options.filter_vcf_no and match.chip_match:
            continue
        else:
            to_be_printed += 1
        if chips != "":
            names += ","
            varscores += ","
            refscores += ","
            chips += ","
        if vargc != "":
            varht += ","
            refht += ","
            vargc += ","
            refgc += ","
        names += match.name
        varscores += str(round(match.var_score, 4))
        refscores += str(round(match.ref_score, 4))
        # Sequence environment data
        if match.var_gc is not None:
            varht += sublist_str(match.var_ht, 4)
            refht += sublist_str(match.ref_ht, 4)
            vargc += str(round(match.var_gc, 4))
            refgc += str(round(match.ref_gc, 4))
        if match.chip_match:
            chips += "Y"
        else:
            chips += "N"

    # If there are no matches, print the line unchanged or filter it out (return
    # without printing)
    if to_be_printed == 0:
        if (not options.filter_vcf_motif and not options.filter_vcf_chip
            and not options.filter_vcf_no):
            print(line, file=output_fileHandle)    # QQQ: ide complains about syntax, but test ran in interactive mode; is this line bad?
                                                   # YYY: Seems fine to me? What's the syntax error? Indentation is a little weird
                                                   #     looking due to the multiline if statement, but should work fine.
        return

    outline = ""
    idx = 0

    for col in columns:
        if outline != "":
            outline += "\t"
        if idx == 7:
            outline += "MOTIFN=" + names + ";MOTIFV=" + varscores + ";MOTIFR=" + refscores
            if refgc != "":
                outline += ";MOTIFVH=" + varht + ";MOTIFRH=" + refht
                outline += ";MOTIFVG=" + vargc + ";MOTIFRG=" + refgc
            if (options.chip_present):
                outline += ";MOTIFC=" + chips
            if col != '.':
                outline += ";" + col
        else:
            outline += col
        idx += 1

    if idx < 7:
        print("**Error** VCF formatted incorrectly. Less than 8 columns found:\n" + line)
        # Add output at the end anyway
        outline += "\tMOTIFN=" + names + ";MOTIFV=" + varscores + ";MOTIFR=" + refscores
        if (options.chip_present):
            outline += ";MOTIFC=" + chips

    print(outline, file=output_fileHandle)

    return


def sublist_str(sublist, sig_figs):
    # Converts a float list within an info field list into a string
    output = ""
    for float_item in sublist:
        if output != "":
            output += "/"
        # debug output += str(float_item)
        output += str(round(float_item, sig_figs))
    if output != "":
        return "(" + output + ")"
    return ""


def chr_less(chr_left, chr_right, sorted_lex):
    """
    Returns true if the left chromosome comes before the right or is the same.

    Args:
        chr_left = (string) the left chromsome being compared
        chr_right = (string) the right chromosome being compared
        sorted_lex = true if file sorted lexicographically (by string values)

    Returns: Whether or not left chromosome comes before the right (boolean).
    """

    # True if the file is sorted lexicographically
    # i.e. chr1 < chr11 < chr2 < chrX < chrY
    if sorted_lex:
        return (chr_left < chr_right)

    # False if the file is sorted numerically
    # i.e. chr1 < chr2 < chr11 < chrX < chrY
    else:
        left = chr_left[3:]    # assumes chromosome name is chr<other bits>
        right = chr_right[3:]
        try:
            l_num = int(left)
            try:
                r_num = int(right)
                return l_num < r_num
            # Right chromosome is a string (chrX, chrY, chrTest, etc)
            except:
                # Left is a number and right is a string
                # Numbers are sorted before strings (chr1 < chrX)
                return True
        # Left number is a string if get to ValueError exception
        except ValueError:
            try:
                r_num = int(left)
                # Left is a string and right is a number
                # Numbers are sorted before strings (chrX !< chr1)
                return False
            # Both are strings, sort lexicographically
            except ValueError:
                return chr_left < chr_right


def match_peaks(chrom, pos, peaks, chip_fh, matches, output_fileHandle, sorted_lex, filter_bed):
    """
    Returns an array of peaks that match the current chromosome and position.
    Updates the output_fileHandle if not None.

    Args:
        p_chr = (string) previous chromosome. Needed to know if a new chromosome
            is being entered.
        chr = (string) chromosome. Chromosomes 1-22, X, and Y are expected.
        pos = current position of the variant.
        peaks = buffer of peaks. They should all be upstream of or overlapping
            variant at chr and pos. Peaks is an array of tuples of the form:
            (chr, start pos, end pos, array of ChIP tfs, motif match array)
            tfs = transcription factors
        chip_fh = input ChIP bed file handle to be read from.
            This must have the start positions in order
            within each chromosome and must be grouped by chromosome.
        matches = list of motif matches as tuples of the form:
            (name, variant score, reference score, ChIP match)
        output_fileHandle = ChIP output bed file to be printed to.
        sorted_lex = True if sorted lexigraphically (by character strings)
        filter_bed = if true exclude non matched motif(s) in ChiP output

    Returns: Peak buffer tuple of the form ( overlapping peak array, next peak )
        Array of peaks that overlap the current chromosome and position
        Next peak (because you must over-read to make sure you don't miss any)
        match_peaks also updates output_fileHandle.
    """

    if chip_fh is None:
        return (peaks, matches)

    # Get rid of peaks that are upstream (left of) of the current chromosome
    idx = 0
    while idx < len(peaks) and chr_less(peaks[idx][0], chrom, sorted_lex):
        # If the chromosome doesn't match, output the line and keep searching
        print_peak(peaks[idx], output_fileHandle, filter_bed)
        idx += 1

    # peak at idx will be included and the rest will be removed
    peaks = peaks[idx:]

    # If previous peaks were not from correct chromosome, get there
    if (len(peaks) == 0):
        new_peak = get_peak_at(chrom, pos, chip_fh, output_fileHandle, sorted_lex, filter_bed)
        # If end of file is reached
        if (new_peak is None):
            return ([], matches)
        else:
            peaks.append(new_peak)

    idx = 0

    # Read through bed file
    while True:

        # If more peaks are needed
        if idx == len(peaks):
            n_peak = get_next_peak(chip_fh)

            # If end of bed file is reached, then just return current list
            if n_peak is None:
                return (peaks, matches)

            peaks.append(n_peak)

        # Current peak (chromosome, start pos, end pos,
        #                transcription factor array, matrix match array)
        #    when get_next_peak defines peak, matrix match array is undefined
        (pchr, psta, pend, ptfs, pmms) = peaks[idx]

        # If next chromosome is reached in bed file [QQQ: should this occur before append?]
        if pchr != chrom:
            break

        if psta <= pos:
            if pend >= pos:
                motif_idx = 0
                for motif_idx in range(len(matches)):
                    pmms.append(matches[motif_idx])    # defines pmms for peaks
                    for trans_factor in ptfs:
                        # If the transcription factor (chip peak) name is the same as
                        # the matched motif name, note that there is a chip match
                        if trans_factor == matches[motif_idx].name:
                            # Motif match is verified by ChIP data
                            matches[motif_idx].chip_match = True
                # Save with new value for pmms
                peaks[idx] = (pchr, psta, pend, ptfs, pmms)
            # Otherwise both are before pos, so remove that peak and continue
            # This should only ever happen when idx is 0... but still
            else:
                print_peak(peaks[idx], output_fileHandle, filter_bed)
                peaks = peaks[0:idx] + peaks[idx + 1:]
                idx -= 1
        # Otherwise peak start is after the variant position, so stop
        else:
            break
        idx += 1

    return (peaks, matches)


def get_next_peak(opened_file):
    """
    Reads in the next line of the bed and returns the next peak's information

    Args: opened_file = an already open input .bed file handle

    Returns: a tuple with the following information (in order) or None
        chromosome number as a string e.g. "chr1",
        position of the start of the peak
        position of the end of the peak
        array list containing transcription factors which bind in that area
        array list containing motif matches (empty)
    """

    # print("Entering get_next_peak")
    # sys.stdout.flush()

    line = opened_file.readline().strip()

    # print("Got: <"+line+">")
    sys.stdout.flush()

    # input file is empty
    if line == "":
        return None

    line_list = line.split('\t')

    chrom = line_list[0]
    start = int(line_list[1])
    end = int(line_list[2])
    tf_array = line_list[3].split(';')

    return (chrom, start, end, tf_array, [])


def get_peak_at(chrom, pos, chip_fh, out_fh, sorted_lex, filter_bed):
    """
    Get the first peak where the end point of the peak is past the input pos.
    Requires that chip_fh is sorted the same way is the vcf input file.
        This file will print all intermediate peaks to out_fh if not None.

    Args:
        chr = the new chromosome to get to
        pos = the new position to get to
        chip_fh = an already open input .bed file handle
        out_fh = an already open output file (to be printed to). May be None.
        sorted_lex = True if sorted lexigraphically (by character strings)
        filter_bed = if true exclude non matched motif(s) in ChiP output

    Returns:
        The first peak from chip_fh that is the same chromosome as chr
        where the end point of the peak is past the input position (pos).
        If the end of file is reached, None is returned.
    """

    # print("*Entering get_peak_at")
    # sys.stdout.flush()

    # Skip ahead until the correct chromosome
    peak = get_next_peak(chip_fh)

    while peak is not None:
        (p_chr, p_sta, p_end, p_tfa, p_mm) = peak
        # If the chromosome is correct, skip ahead until the correct position
        if p_chr == chrom:
            # We have passed the position at the chromosome
            if p_end >= pos:
                # print("get_peak_at returns: "+p_chr+":"+str(p_sta)+"-"+str(p_end))
                # sys.stdout.flush()
                return peak
            else:
                print_peak(peak, out_fh, filter_bed)

        # If the chromosome is too low and there is an outfile, print
        elif chr_less(p_chr, chrom, sorted_lex):
            print_peak(peak, out_fh, filter_bed)
        # If we have passed the chromosome
        else:
            return peak

        peak = get_next_peak(chip_fh)
    # print("*get_peak_at returns None")
    # sys.stdout.flush()
    return None


def print_peak(peak, fileHandle, filter_bed):
    """
    Prints the peak to the given file handle (or exits if no file is given)

    Args:
        peak = ChIP peak of the form:
            (chr, start, stop, chip tf array, motif match tf array)
            chip tf array is an array of tf names
            motif match tf array is an array of (motif name, vscore, rscore, strand)
        fileHandle = the file to print to. This file should already be opened.
            If this file is None, nothing will happen.
        filter_bed = if true exclude lines in chip (bed) output file if they
                     don't match a motif

    Returns: Nothing
    """
    if fileHandle is None:
        return

    (chrom, start, end, c_array, mm_array) = peak

    # If there are no motif matches and filtering is on, do not print this peak
    if filter_bed and len(mm_array) == 0:
        return

    line = chrom + '\t' + str(start) + '\t' + str(end) + '\t'

    # Generate string of chip transcription factors from the peak
    chip_string = ""
    for tf_name in c_array:
        if chip_string != "":
            chip_string += ";"
        chip_string += tf_name

    # Generate string of motif matches that overlap the peak
    motif_string = ""
    for match in mm_array:
        if motif_string != "":
            motif_string += ";"
        motif_string += match.name + "," + str(round(match.var_score, 4)) + ","
        motif_string += str(round(match.ref_score, 4))

    line += chip_string + "\t" + motif_string
    print(line, file=fileHandle)

    return

"""
    ---- START OF MAIN ----

    Requires that vcf file have variants sorted by position within chromosomes
"""


# Create arguments and options
parser.add_argument("-i", "--input", dest="input_file", required=True)
parser.add_argument("-r", "--ref", dest="file_reference", required=True)
parser.add_argument("-m", "--motif", dest="file_motif", required=True)
parser.add_argument("-o", "--output", dest="output_file", required=True)
parser.add_argument("-ci", "--chip", dest="chip_file",
                    required=False, default=None)
parser.add_argument("-co", "--chipout", dest="chip_out_file",
                    required=False, default=None)
parser.add_argument("-bp", "--baseline", dest="baseline_file",
                    required=False, default=None)
parser.add_argument("-th", "--threshold", dest="threshold",
                    required=False, default=None)
parser.add_argument("-pc", "--pseudocounts", dest="pseudocounts",
                    required=False, default=0.1)
parser.add_argument("-ws", "--wing_size", dest="wing_size",
                    required=False, default=50)
parser.add_argument("-fm", "--filter_motif", action="count", required=False)
parser.add_argument("-fc", "--filter_chip", action="count", required=False)
parser.add_argument("-fn", "--filter_novel", action="count", required=False)
parser.add_argument("-fp", "--filter_co", action="count", required=False)
parser.add_argument("-sk", "--kary_sort", action="count", required=False)
parser.add_argument("-mv", "--multi_variant", dest="multi_var",
    required=False, default=-1)    # use -1 and correct to wing_size below
parser.add_argument("-rf", "--force_ref_match", dest="force_ref_match",
    required=False, default=False)

args = parser.parse_args()

# Easier to use argument variables (YYY: why track duplicate variable names?)
# YYY: Was probably just copy and pasted, tbh. No specific reason they were duped.
file_input = args.input_file
file_reference_genome = args.file_reference
file_motif = args.file_motif
file_output = args.output_file
file_baseline_prob = args.baseline_file    # defaults to None
file_chip = args.chip_file
file_output_chip = args.chip_out_file
filter_bed = False
pc = float(args.pseudocounts)
ws = int(args.wing_size)
if (args.multi_var is None):
    multivar_computation_flag = False
else:
    if (int(args.multi_var) == -1):
        multivar_distance = ws
        multivar_computation_flag = False
    else:
        multivar_distance = int(args.multi_var)
        multivar_computation_flag = True

# Options list. Easier to pass in methods or use in code updates.
options = Options_list()
options.filter_vcf_motif = (args.filter_motif is not None)

# Are input files chr sorted lexicographically (or by karyotype order)?
# Input vcf file and chip bed file must be sorted the same way
# -sk sets this to True (means sorted lexicographically)
sorted_lex = (args.kary_sort is None)

print("Run started at:" + timeString())

# Output so user can double check options
print(("Input file: {}\nReference file: {}\nMotif file: {}\n" +
       "Output file: {}\nOptional arguments:\n    Pseudocounts value = {}").
      format(
    file_input, file_reference_genome, file_motif,
    file_output, pc
))

# Optional arguments to print
opt_args = ""

if args.threshold is not None:
    th = float(args.threshold)
    opt_args += "    Defined match threshold = " + str(th) + "\n"
else:
    th = 0.0
    opt_args += "    Default match threshold = " + str(th) + "\n"

fileHan_chip = None
fileHan_out_chip = None

if (file_chip is not None):
    options.chip_present = True

    # -- Process chip related options
    # Should lines in the chip (bed) output file be excluded
    # if they don't match a motif?
    # -fp sets this to True
    filter_bed = (args.filter_co is not None)

    options.filter_vcf_chip = (args.filter_chip is not None)
    options.filter_vcf_no = (args.filter_novel is not None)

    opt_args += "    ChIP file: " + file_chip + "\n"
    fileHan_chip = open(file_chip)
    if (file_output_chip is not None):
        opt_args += "    ChIP output file: " + file_output_chip + "\n"
        fileHan_out_chip = open(file_output_chip, "w")
        if (filter_bed):
            opt_args += "    Filter output ChIP bed for motif matches? Yes\n"
        if (options.filter_vcf_chip):
            opt_args += "    Filter output vcf for ChIP peak overlap? Yes\n"
        if (options.filter_vcf_no):
            opt_args += "    Filter output vcf for no ChIP peak overlap? Yes\n"
        if (options.filter_vcf_chip and options.filter_vcf_no):
            # YYY: default to one? or should one override the other?
            opt_args += "Err: Cannot have -fn and -fc (two prev options)."
            opt_args += "    Both will be ignored.\n"
            options.filter_vcf_chip = False
            options.filter_vcf_no = False
elif (file_output_chip is not None):
    opt_args += "No ChIP file given, so no ChIP output file will be created\n"

if (file_baseline_prob is not None):
    opt_args += "    Baseline probabilities file: " + file_baseline_prob + "\n"
if (options.filter_vcf_motif):
    opt_args += "    Filter output vcf for motif matches? Yes\n"

if (not sorted_lex):
    opt_args += "    Input vcf and Input ChIP bed are sorted by karyotype\n"
# XXX: add multivariant option output and force reference
print(opt_args)    # debug outputs

# Main -- Post Arguments

if (file_baseline_prob is None):
    bp = [0.25, 0.25, 0.25, 0.25]
else:
    print("Reading in baseline probabilities (@" + timeString() + "):")
    bp = motif.get_baseline_probs(file_baseline_prob)
    print(bp)
sys.stdout.flush()

# Grab motif list from motif file.
print("Creating motif list from " + format(file_motif) +
      ", with pc=" + format(pc) +
      ", threshold=" + format(th) +
      ", bp=" + format(bp))
sys.stdout.flush()
#(motif_set, max_motif_l) = get_motifs(file_motif, pc, th, bp)    # old call
motif_set = motif.get_motifs(file_motif, pc, th, bp)
    # motif_set.max_positions              replaces max_motif_l
    # motif_set.motifs[index].positions    more 'valid'
# debug print("Maximum motif length is "+str(motif_set.max_positions)+".")

# need to grab sequence as full variant but only pull in
# length around wing that is necessary!
#     XXX:QQQ: what does it mean to use wing_l
#        as the max of longest motif and argument ws?
#        if either > individual motif length --> only uses motif length
#        therefore wing_l can not be greater than motif length
#        if wing_l less than max motif length then larger motifs are
#        ignored during processing? what happens? QQQ
# Wing length for individual variants is calculated dynamically so that it is
# always one less than the length of the maximum motif so motif is only matched
# against places with overlap with the variant
# QQQ: should wing_length be the number of elements to grab or the maximum index?
#        ie motif_set.max_positions or motif_set.max_positions - 1
wing_l = min( motif_set.max_positions, ws ) #max(motif_set.max_positions - 1, ws)
print("Proceeding with wing_size: " + format(wing_l) + " vs defined ws(" + format(ws) + ")" )

# Open output file.
fileHan_output = open(file_output, "w")
# XXX: should check for output open here

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
print("Creating index from reference sequence for efficient searching..." +
      timeString())
print("This will be slow the first time (takes about 20 seconds on i7.)")
sys.stdout.flush()
fa_ind = Fasta(file_reference_genome)    # XXX: need to check, if present skip
print("Completed fasta index @ " + timeString())

print("Importing variants(" + timeString() + "). This may take a while.\n")
sys.stdout.flush()

# Open and Read VCF file: populates SequenceArray, assumes set fits in memory
variant_set = sequence.SequenceArray()
with open(file_input) as vcf_handle:

    line = vcf_handle.readline()

    # building output vcf info line
    info_needed = True
    info = "##INFO=<ID=MOTIFN,Number=.,Type=String,Description="
    info += "\"Matched motif names\">"
    info += "\n##INFO=<ID=MOTIFV,Number=.,Type=Float,Description="
    info += "\"Variant match scores\">"
    info += "\n##INFO=<ID=MOTIFR,Number=.,Type=Float,Description="
    info += "\"Reference match scores\">"

    info += "\n##INFO=<ID=MOTIFVH,Number=.,Type=String,Description="
    info += "\"Variant homotypic match scores\">"
    info += "\n##INFO=<ID=MOTIFRH,Number=.,Type=String,Description="
    info += "\"Reference homotypic match scores\">"
    info += "\n##INFO=<ID=MOTIFVG,Number=.,Type=String,Description="
    info += "\"Variant environment GC content\">"
    info += "\n##INFO=<ID=MOTIFRG,Number=.,Type=String,Description="
    info += "\"Reference environment GC content\">"

    if (fileHan_chip is not None):
        info += "\n##INFO=<ID=MOTIFC,Number=.,Type=Character,Description="
        info += "\"Motif validated by ChIP (Y/N)\">"

    print("\tfinished header read " + timeString())

    # Skip info lines
    while line.startswith("##"):
        # Print new info lines at the top of the ##INFO section
        if info_needed and line.startswith("##INFO"):
            print(info, file=fileHan_output)
            info_needed = False
        print(line, file=fileHan_output, end="")
        line = vcf_handle.readline()

    # Create appropriate header.    Presumably reads the first # line
    header = line.strip()
    print(header, file=fileHan_output)
    # QQQ|YYY: push header line to global variable for samples, also modify
    #    info line as method above to all extension to calls and reduce file
    #    read?
    # YYY: Not entirely sure what your question is asking here, but using a
    #     method to add the INFO fields metadata in the header would be fine
    #     and probably look a lot cleaner.
    # CCC-WK: note has to do with whether it used again or not, best way to handle

    # Process each variant
    eof_counter = 0
    while True:
        # Reads in the next line of the vcf
        # adds the next variant's information to the element array object

        line = vcf_handle.readline()

        if line is None:    # stop infinite loop
            print('debug: stop file read with count' + format(eof_counter))
            break

        # skip empty and information lines
        if line.startswith("##"):
            continue
        if line == "":
            # handle while True infinite loop try to find End of File
            if (eof_counter > 1):
                break
            else:
                eof_counter += 1

        line = line.strip()

        # skip empty lines (after removing arbitrary whitespace characters)
        if line == "":
            continue

        # -- process valid variant lines
        line_list = line.split("\t")
        # variant creation: add_seq_defined(chromosome, position, reference_seq, variant_seq):
        # variant_set.add_seq_defined(line_list[0], int(line_list[1]), line_list[3], line_list[4])
        new_sequence_element = sequence.SequenceElement()
        new_sequence_element.assign(line_list[0], int(line_list[1]), line_list[3], line_list[4])
        # grab samples for variant
        new_sequence_element.assign_samples(line_list[9:])
        # push full line (memory hog, but allows multivar computation
        #    w/o significant code manipulation to track which line is current
        #    already processed, etc.)
        new_sequence_element.vcf_line = line
        variant_set.add_seq(new_sequence_element)
        print("\tadd element " + new_sequence_element.name + "::" + timeString())


print("Finished importing variants(" + timeString() + ")\n")

if multivar_computation_flag:
    print("Start variant merge (" + timeString() + ").\n")
    variant_set.multivariant_list_build(multivar_distance)
    print("Finished variant merge(" + timeString() + ").\n")
    # XXX: sort variant_set by the chromosome name? karotype? use chr_less
    #     or rather just insert multivariant items into middle of list
    #     before next chromosome elements
    # YYY: Typically, we'd expect lexicographic order for input, it's the standard
    #     order for sorting this data typically. I don't think we need to build in
    #     a sort, just tell users to sort their VCFs/bed files beforehand.
    # CCC-WK: not user input. the sort has to do with creating multivariants from the
    #     set of individual inputs. The program is doing work for the user, by
    #     combining variants into super-variants, see email discussion.
    #     so there needs to be a decision on where in the order to put
    #     the super-variants for processing. also how to annotate in output

print("Analyzing variants(" + timeString() + "). This may take a while.\n")

# ---------------------------------------------------------------------------

# Queue of ChIP peaks that overlap the current variant
# Will contain peaks as a tuple of the form
# (chr, start, stop, chip tf array, motif match tf array)
# chip tf array is an array of tf names
# motif match tf array is an array of (motif name, vscore, rscore, strand)
peak_buffer = []
chromosome = ""    # processed chromosome
for index in range(variant_set.length()):
    var_element = variant_set.seq[index]

    # QQQ: separate variant and multivariant test sets?
    #    test for and only analyze if index not in variant_set.multivariant
    # YYY: Seems this will likely be a necessity, at least to me. May be a way
    #     around it, but I can't think of a better option off the top of my head.

    # Update previous, next, and current variables
    if (chromosome != var_element.name):
        chromosome = var_element.name
        print("\tStart Analyzing new chromosome: " + chromosome + "...")

    # -- prepare the variant for processing
    # 1. Get reference sequence surrounding the variant from the reference file
    #    get_surrounding_seq(chr, pos, len(ref_bases), wing_l, fa_ind)
    #    self.get_surround_seq(self, wing_length, fasta_object, force_ref_match)
    #
    print("\t\tpulling " + format(wing_l) +
        " base wings from index " + format(fa_ind) +
        " arg: " + format(args.force_ref_match))
    var_element.get_surround_seq(wing_l, fa_ind, args.force_ref_match)
    # 2. compute reverse complement
    var_element.assign_rev_complement()
    # 3. compute int version (XXX: after testing for speed up is confirmed)

    print("Calculating:\n" + var_element.print_str() + "at(" + timeString() + ")\n")

    # -- building the list of MotifMatch objects
    # Calculate motif matches to variant sequence
    ref_seq = var_element.return_full_ref_seq(wing_l)
    var_seq = var_element.return_full_var_seq(wing_l)
    plusmatch = motif_set.motif_match(bp, ref_seq, var_seq, wing_l)
    #     plusmatch returns an list of MotifMatch objects
    # Add local environment data: XXX: new version guessed intent see code
    #    the process_local_env function never returned anything before
    #    the individual iteration output was accumulated in a transient variable
    #    and the output (presuming PLE had any) was not previously assigned
    #
    plusmatch = motif_set.process_local_env(bp, plusmatch, var_element, None, var_seq, ref_seq)

    # Calculate motif matches to reverse complement
    ref_seq_rc = var_element.return_full_ref_seq_reverse_complement(wing_l)
    var_seq_rc = var_element.return_full_var_seq_reverse_complement(wing_l)
    minusmatch = motif_set.motif_match(bp, ref_seq_rc, var_seq_rc, wing_l)
    # Add local environment data
    minusmatch = motif_set.process_local_env(bp, minusmatch, var_element, None, var_seq_rc, ref_seq_rc)

    matches = plusmatch + minusmatch
    print(("\t" + format(var_element.name) + ":" + format(var_element.position) +
        " +match(" + format(len(plusmatch)) +
        ") + -match(" + format(len(minusmatch)) +
        ") = matches(" + format(len(matches)) + ")"))

    # Update ChIP buffer for current position
    # Update matches array with peak overlap data
    #    WARNING: XXX: match_peaks has not been heavily reviewed for validity
    #    note: fileHan_chip is only read by children of match_peaks
    (peak_buffer, matches) = match_peaks(var_element.name, var_element.position,
                                         peak_buffer, fileHan_chip,
                                         matches, fileHan_out_chip,
                                         sorted_lex, filter_bed)

    print(("matches number:" + format(len(matches))))
    """print("match_peaks returned "+str(len(peak_buffer))+" peak(s):")
    for peak in peak_buffer:
        (pchr, psta, pend, ptfs, pmms) = peak
        print(pchr+":"+str(psta)+"-"+str(pend)+" tfs:"+str(len(ptfs))+
            " mms:"+str(len(pmms)))
    print()"""

    # Co-binding transcription factors currently not implemented
    cb_dict = None

    # Create the correct line in VCF format and print to file_output
    update_vcf(var_element.vcf_line, matches, fileHan_output, options)
    sys.stdout.flush()

    # Print remaining peaks
    for peak in peak_buffer:
        print_peak(peak, fileHan_out_chip, filter_bed)


print("Finished analyzing variants(" + timeString() + ").\n")

# Close output files.
fileHan_output.close()
if fileHan_out_chip is not None:
    fileHan_out_chip.close()
if fileHan_chip is not None:
    fileHan_chip.close()
