#!/usr/bin/env python
"""
Calculate motif thresholds for a given motif file and print to a new file.

Usage: thresholds.py -m <motifs.txt> -o <output.txt> [OPTIONS]

Args:
    -m (required) = Input filename of a file containing PWMs.
    -o (required) = Output filename.
    -bp (optional) <baselines.txt> = A file containing a single line with tab delineated values for baseline
        probabilities for A, C, G, T (in order). Probabilities should all be positive and should sum to 1. If not
        provided then all are assumed to be equally likely (all are 0.25).
    -pc (optional) <0.1> = Pseudocounts value to be added to all positions of the motif frequency matrix before
        calculating the probability matrix.
    -th (optional) <0.0> = Default threshold value. This is used if the calculated threshold is lower.
        Ex: default_th = 0.0, biopython calculates threshold needed for a given false positive rate is -1.23, threshold
        printed will be 0.0.
    -fpr (optional) <0.05> = Acceptable false positive rate for defining thresholds for each motif.
    -pe (optional) <4> = Precision exponent used for threshhold calculations. Making this greater than 5 may result in
        extremely slow run times. Using a lower number will result in faster (but potentially innacurate) calculations.
    -ow (optional flag) = OverWrite: If present, thresholds already present in the input file will be replaced in the
        output file.
"""
from Bio import motifs
import sys
import argparse
parser = argparse.ArgumentParser(usage=__doc__)


def get_baseline_probs(baseline_f):
    """
    Read in baseline probabilities from a file that has them listed in the first line.

    Args:
        baseline_f (str): Name of file containing a probability array of the form:
                [ PrA PrC PrT PrG ] where PrA + PrC + PrT + PrG = 1 (and all are positive and non-zero).

    Returns:
        bp_array (list): List of the probabilities for each base: [ PrA, PrC, PrT, PrG ]
    """

    # Default baseline probability numbers (assumes all are equally likely).
    bp_array = [0.25, 0.25, 0.25, 0.25]
    if baseline_f is None:
        return bp_array

    with open(baseline_f) as f:
        try:
            for line in f:
                # remove commas, brackets, and whitespace on far left and right
                line = line.strip().replace(
                    '[', '').replace(']', '').replace(',', '')
                if line != "":
                    line = line.split()
                    for idx in range(4):
                        bp_array[idx] = float(line[idx])
                    return bp_array
        except ValueError:
            print("**ERROR** Baseline probability file incorrectly formatted.",
                  "\nFile should contain only [ PrA PrC PrT PrG ]\n",
                  "Where PrA + PrC + PrT + PrG = 1 (and all are positive and",
                  "non-zero)\n",
                  "Continuing with:")
            return bp_array

        print("**ERROR** Empty baseline probability file found.\n",
              "File should contain only [ PrA PrC PrT PrG ]\n",
              "Where PrA + PrC + PrT + PrG = 1 (and all are positive and",
              "non-zero)\n",
              "Continuing with:")
        return bp_array


def output_motifs(input_f, output, default_th, overwrite, thresholds_list):
    """
    Read in and calculate probability matrices from a frequency matrix file.

    Args:
        input_f (str): Name of file containing frequency matrices for each motif.
        output (str): Name of output file.
        default_th (float): Default threshold value. This is used if the calculated threshold is lower than this value.
            This value may be None.
            Ex: default_th = 0.0, biopython calculates threshold needed for
            a given false positive rate = -1.23, threshold printed will be 0.0.
        overwrite (bool): True if thresholds already in the file should be replaced.
        thresholds_list (list): List of thresholds (floats) calculated by biopython.
    """

    output_f = open(output, "w")
    idx = 0

    with open(input_f) as f:

        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n
        # arrays for C, G, T - each with same format as A
        iden = "No id found"
        name = "No name found"
        # Index for line in the motif matrix.
        i = 0
        # Position frequency matrix (read in then reprinted unchanged)
        pfm = ["", "", "", ""]
        given_thresh = None

        for line in f:

            # First line contains id and name
            if i == 0:
                line = line.strip().split()
                iden = line[0]
                name = line[1]
                # Read in listed threshold
                if (len(line) > 2):
                    given_thresh = float(line[2])

            # Order of position weight matrices is A,C,G,T.
            elif i < 5:
                pfm[i - 1] = line

            # Output motif and continue (there are 2 newlines between motifs)
            else:
                if not overwrite and given_thresh is not None:
                    th = given_thresh
                elif (default_th is not None and
                      default_th > thresholds_list[idx]):
                    th = default_th
                else:
                    th = thresholds_list[idx]
                # print the motif back out to the output file
                out_line = iden + "\t" + name + "\t" + str(th) + "\n"
                out_line += pfm[0] + pfm[1] + pfm[2] + pfm[3]
                print(out_line, file=output_f)
                idx += 1
                i = -1
                given_thresh = None
                pfm = ["", "", "", ""]
            i += 1

    if i >= 5:
        if not overwrite and given_thresh is not None:
            th = given_thresh
        elif default_th is not None and default_th > thresholds_list[idx]:
            th = default_th
        else:
            th = thresholds_list[idx]
        # print the motif back out to the output file
        out_line = name + "\t" + iden + "\t" + str(th) + "\n"
        out_line += pfm[0] + pfm[1] + pfm[2] + pfm[3]
        print(out_line, file=output_f)

    output_f.close()
    return

# TODO - Main statement, etc.
parser.add_argument("-m", "--motif", dest="motif_file", required=True)
parser.add_argument("-o", "--outfile", dest="motif_outfile", required=True)
parser.add_argument("-bp", "--baseline", dest="baseline_file", required=False, default=None)
parser.add_argument("-pc", "--pseudocounts", dest="pseudocounts", required=False, default=0.1, type=float)
parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=None)
parser.add_argument("-fpr", "--falsepos", dest="false_pos_rate", required=False, default=0.05, type=float)
parser.add_argument("-pe", "--precision", dest="precision_exp", required=False, default=4, type=int)
parser.add_argument("-ow", "--overwrite", action="store_true", required=False)

args = parser.parse_args()

if args.baseline_file is not None:
    bp = [0.25, 0.25, 0.25, 0.25]
else:
    bp = get_baseline_probs(args.baseline_file)
pc = args.pseudocounts
if args.threshold is not None:
    d_th = float(args.threshold)
else:
    d_th = None
ow = args.overwrite is not None
fpr = float(args.false_pos_rate)
pe = int(args.precision_exp)
if pe > 5:
    print("Warning: high precision exponent (-pe) may cause drastic slowing " + "or memory errors")
if pe <= 0:
    pe = 1
    print("Precision exponent (-pe) too low, set to " + str(pe))

thresholds = []
background = {'A': bp[0], 'C': bp[1], 'T': bp[2], 'G': bp[3]}
print("Baseline nucleotide frequencies:\n" + str(background))


print("Calculating thresholds. This could take a while.")
sys.stdout.flush()
idx = 0
exponent = 1

# Calculate thresholds using biopython
fh = open(args.motif_file)
for m in motifs.parse(fh, "jaspar"):
    pwm = m.counts.normalize(pseudocounts=pc)
    pssm = pwm.log_odds(background)
    # Precision argument of 4 was recommended by biopython's documentation
    distribution = pssm.distribution(background=background, precision=10**pe)
    m_th = distribution.threshold_fpr(fpr)
    thresholds.append(m_th)
    # print progress
    idx += 1
    if (idx >= 10**exponent):
        print(str(idx) + " thresholds calculated...")
        exponent += 1
        sys.stdout.flush()

print("Total motifs read: " + str(len(thresholds)))

print("Outputing thresholds")
output_motifs(args.motif_file, args.motif_outfile, d_th, ow, thresholds)

print("Done.")
