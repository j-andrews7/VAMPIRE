#!/usr/bin/env python
"""
For a given motif annotated vcf file (already run through motifs.py), remove all motif matches for 
TFs that are not expressed in at least one sample above the specified threshold.
    
Usage: activity.py -i <input.vcf> -e <expression.bed> -o <output.vcf> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process. 
    -o (required) <output.vcf>: Name of output file to be created.
    -e (required) <expression.bed>: An expression 'bed' file.
    -th (optional) <5>: TFs are considered expressed if they are above this threshold. 
"""
# TODO - Write doc string, usage statement, etc. Clean up old comments, debug statements, etc.

import sys
import argparse
from statistics import mean, stdev

class Position:
    """
    Use to represent and handle genomic ranges more easily.

    Args:
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
    """

    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start = start_pos
        self.end = end_pos


    def overlaps(self, pos_b):
        """
        Return whether self overlaps Position pos_b.

        Args: pos_b = (Position object) another position

        Returns: 
            (bool): True if self overlaps with Position pos_b. False if not.
        """
        if pos_b == None:
            return False

        if self.chrom != pos_b.chr:
            return False

        start_max = max(self.start, pos_b.start)
        end_min = min(self.end,   pos_b.end)

        return start_max <= end_min


def parse_activity_file(activity_file):

    with open(activity_file) as f:
        header = f.readline.strip()
        # TODO - Parse header to get sample names and indices.
        # TODO - Determine how to return enhancer info. dict or whatever.

        for line in f:
            line_list = line.strip().split('\t')

            # Create position object for enhancer.
            chrom = line_list[0]
            start_pos = int(line_list[1])
            end_pos = int(line_list[2])
            pos = Position(chrom, start_pos, end_pos)

            iden = line_list[3]

            # Get activity scores for samples
            samples_act = [float(x) for x in line_list[4:]]

            return (pos, iden, samples_act)


def output_activity(open_file, enh_pos, enh_id, matches, ref_l, var_l, options):

    # If filtering and there are no matches, print nothing
    if options.filter_bed and len(matches) == 0:
        return

    # This is set to 0 by default, so everything will be output
    if len(matches) < options.filter_bed_num:
        return

    line = enh_pos.chrom + '\t' + str(enh_pos.start) + '\t' + str(enh_pos.end) + '\t' + enh_id
    # print number of samples affected
    if len(matches) > 0:
        line += "\tsig_diff=" + str(len(matches)) + ",have_var=" + str(var_l)
        line += ",ref_group=" + str(ref_l)
    else:
        line += "\t."

    for match in matches:
        (name, pos, zscore, orig_line) = match

        # info line is 7th column
        line_list = orig_line.split("\t")
        if len(line_list) > 7:
            out_str = ""
            for field in line_list[7].split(';'):
                if field.startswith("MOTIF"):
                    out_str += field + ";"
            # Cut off last ';'
            out_str = out_str[:-1]

        line += '\tsample=' + name + ",var=" + str(pos) + ",z=" + str(round(zscore, 3)) + "--" + out_str

    print(line, file=open_file)


####-PARSER-####
parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument("-a", "--activity", dest="activity_file", required=True)
parser.add_argument("-o", "--output", dest="output_file", required=True)
parser.add_argument("-i", "--input", dest="input_file", required=True)
parser.add_argument("-th", "--threshold", dest="threshold",
                    required=False, default=0)
parser.add_argument("-fan", "--filter_a_n", dest="filter_a_n",
                    required=False, default=0)         
        # Only print activity if it affects more samples than this number
        # default is -1 so a region with any number of samples affected
        # (having z-scores above threshold) will be output.
parser.add_argument("-eo", "--errout", dest="err_out",
                    required=False, default=None)
parser.add_argument("-fa", "--filter_a", action="count", required=False)          
        # Should lines in the activity (bed) output file be excluded
        # if they don't match a motif?
        # -fa sets this to True
parser.add_argument("-fv", "--filter_vcfs", action="count", required=False)         
        # Should lines in the vcf output files be excluded
        # if they don't have activity?
        # -fv tag sets this to True


args = parser.parse_args()

# Easier to use argument variables
act_file = args.activity_file
out_file = args.output_file
err_file = args.err_out
inp_file = args.input_file
th = float(args.threshold)

options = Options_list()
options.filter_bed_num = int(args.filter_a_n)
options.filter_bed = (args.filter_a != None)
options.filter_vcf = (args.filter_vcfs != None)


print("Done")
