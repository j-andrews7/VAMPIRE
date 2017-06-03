#!/usr/bin/env python
"""
Calculate motif thresholds for a given motif file and print to a new file.

Usage: thresholds.py -m <motifs.txt> -o <output.txt> [OPTIONS]

Args:
    -m (required) = Input filename of a file containing PWMs.
    -o (required) = Output filename.
    -bp (optional) <baselines.txt> = A file containing a single line with
        tab delineated values for baseline probabilities for A, C, G, T (in order).
        Probabilities should all be positive and should sum to 1. If not
        provided then all are assumed to be equally likely (all are 0.25).
    -pc (optional) <0.1> = Pseudocounts value to be added to all positions of
        the motif frequency matrix before calculating the probability matrix.
    -th (optional) <None> = Default threshold value. This is used if the
        calculated threshold is lower. Good unknown value is zero.
        Ex: for default_th = 0.0, if biopython calculates threshold needed for
        a given false positive rate as -1.23, threshold printed will be 0.0.
    -fpr (optional) <0.05> = Acceptable false positive rate for defining thresholds for each motif.
    -pe (optional) <4> = Integer precision exponent used for threshhold calculations.
        Making this greater than 5 may result in extremely slow run times.
        Using a lower number will result in faster (but potentially innacurate) calculations.
        Allowed value range 1 to Inf (technically allowed, not-advised).
    -ow (optional flag) = OverWrite: If present, thresholds already present in
        the input file will be replaced in the output file.
"""

import sys
import argparse

from Bio import motifs

import motif
from utils import timeString

parser = argparse.ArgumentParser(usage=__doc__)


def main(motif_file, motif_outfile, d_th, pc, bp, ow, fpr, pe):
    thresholds = []

    # -- pre-processing passed variables
    if bp is None:
        bp = [0.25, 0.25, 0.25, 0.25]
    else:
        bp = motif.get_baseline_probs(bp)
    if d_th is not None:
        d_th = float(d_th)
    else:
        d_th = None
    fpr = float(fpr)
    pe = int(pe)
    if pe > 5:
        print(("Warning: high precision exponent (-pe=" +
               str(pe) + ") may cause drastic slowing or memory errors"))
    if pe <= 0:
        pe = 1
        print(("Precision exponent (-pe) too low, set to " + str(pe)))

    # -- start function in earnest
    background = {'A': bp[0], 'C': bp[1], 'T': bp[2], 'G': bp[3]}
    print(("Baseline nucleotide frequencies:\n\t" + str(background)))

    print(("Calculating thresholds (" + timeString() + "). This could take a while."))
    sys.stdout.flush()
    idx = 0
    print_exponent = 1

    # Calculate thresholds using biopython
    fh = open(motif_file)
    for m in motifs.parse(fh, "jaspar"):
        pwm = m.counts.normalize(pseudocounts=pc)    # creates dictionary like representation
        pssm = pwm.log_odds(background)              # converts to log_odds vs background
        # Precision argument of 4 was recommended by biopython's documentation (slow step)
        # YYY - JA - 05/05/2017 - This could be sped up by using TFM-PVALUE's C++ functions for
        # determining thresholds. May consider implementing, calculates precise p-values with no
        # errors much more quickly and can also generate p-values from a score.
        distribution = pssm.distribution(background=background, precision=10 ** pe)
        m_th = distribution.threshold_fpr(fpr)
        thresholds.append(m_th)
        # print progress
        idx += 1
        if (idx >= 10 ** print_exponent):
            print((str(idx) + " thresholds calculated... at " + timeString()))
            print_exponent += 1
            sys.stdout.flush()

    print(("Total motifs read: " + str(len(thresholds))))

    print("Outputing thresholds")
    motif.get_put_motifs(motif_file, motif_outfile, d_th, ow, thresholds)

    print(("Done (" + timeString() + ")"))

if __name__ == '__main__':

    parser.add_argument("-m", "--motif", dest="motif_file", required=True)
    parser.add_argument("-o", "--outfile", dest="motif_outfile", required=True)
    parser.add_argument("-bp", "--baseline", dest="baseline_file", required=False, default=None)
    parser.add_argument("-pc", "--pseudocounts", dest="pseudocounts",
                        required=False, default=0.1, type=float)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False)
    parser.add_argument("-fpr", "--falsepos", dest="false_pos_rate",
                        required=False, default=0.05, type=float)
    parser.add_argument("-pe", "--precision", dest="precision_exp",
                        required=False, default=4, type=int)
    parser.add_argument("-ow", "--overwrite", action="store_true", required=False)

    args = parser.parse_args()

    motif_file = args.motif_file
    motif_outfile = args.motif_outfile
    d_th = args.threshold
    pc = args.pseudocounts
    bp = args.baseline_file
    ow = args.overwrite
    fpr = args.false_pos_rate
    pe = args.precision_exp

    main(motif_file, motif_outfile, d_th, pc, bp, ow, fpr, pe)
