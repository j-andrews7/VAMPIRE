#!/usr/bin/env python
"""
Calculate motif thresholds for a given motif file and print to a new file.

Usage: thresholds.py -m <motifs.txt> -o <output.txt> [OPTIONS]

Args:
    -m (str):
        Input filename of a file containing PWMs.
    -o (str):
        Output filename.
    -a (float, optional) <0.25>:
        Background probability for A nucleotides. If not
        provided then all are assumed to be equally likely (all are 0.25).
    -t (float, optional) <0.25>:
        Background probability for T nucleotides.
    -c (float, optional) <0.25>:
        Background probability for C nucleotides.
    -g (float, optional) <0.25>:
        Background probability for G nucleotides.
    -p (int, optinal) <1>:
        Processor cores to utilize. Will decrease computation time linearly.
    -pc (float, optional) <0.1>:
        Pseudocounts value to be added to all positions of the motif frequency matrix before calculating
        the probability matrix.
    -pv (float, optional) <0.00001>:
        P-value to be used for defining thresholds for each motifs. I don't recommend changing this.
    -ow (flag, optional):
        OverWrite: If present, thresholds already present in the input file will be replaced in the output file.
"""

import sys
import argparse
import gc
from multiprocessing.dummy import Pool as ThreadPool
import itertools
from timeit import default_timer as timer

from Bio import motifs

import readline  # Often necessary for rpy2 to work properly.
import rpy2.robjects as R
from rpy2.robjects.packages import importr

from motif import get_put_motifs
from utils import timeString, check_r_install


def find_thresh(r_func, matrix, pval, r_bg):
    """
    Utilize the given R function to calculate the detection thresholds for each pwm matrix given
    the nucleotide background frequencies and desired p-value.

    Args:
        r_func (R function):
            The R function that will be used to actually calculate the threshold using a matrix,
            p-value, and background nucleotide frequencies.
        matrix (R matrix):
            The PWM for which a detection threshold will be found.
        pval (float):
            P-value to which threshold will be calculated.
        r_bg (R FloatVector):
            Vector containing background nucleotide frequencies [A, C, G, T]
    """
    start = timer()
    thresh = r_func(matrix[1], pval, r_bg)
    end = timer()

    print(matrix[0] + ": " + str(end - start))

    # Attempt to free memory by running R and python garbage collectors. This was an issue.
    R.r("gc()")
    gc.collect()

    return thresh[0]


def main(motif_file, motif_outfile, pc, bp, ow, pv, p):
    matrices = []
    background = {'A': bp[0], 'C': bp[1], 'G': bp[2], 'T': bp[3]}
    r_background = R.FloatVector((bp[0], bp[1], bp[2], bp[3]))
    print(("Baseline nucleotide frequencies:\n\t" + str(background)))

    # Calculate thresholds using biopython.
    print(("Reading in motifs."))
    fh = open(motif_file)
    for m in motifs.parse(fh, "jaspar"):
        pfm = m.counts.normalize(pseudocounts=pc)    # Create frequency matrix.
        pwm = pfm.log_odds(background)              # Calculate to log likelihoods vs background.

        # Create R matrix from motif pwm.
        mat = R.r.matrix(R.FloatVector(pwm[0] + pwm[1] + pwm[2] + pwm[3]), nrow=4)
        matrices.append((m.name, mat))

    fh.close()

    # Define function to calculate thresholds with TFMPvalue R package.
    get_thresh = R.r('''
    function(mat, pvalue, bg) {
        gc()
        rownames(mat, do.NULL = TRUE, prefix = "row")
        rownames(mat) <- c("A","C","G","T")
        bg <- c(A=bg[1], C=bg[2], G=bg[3], T=bg[4])
        TFMpv2sc(mat, pvalue, bg, type="PWM")
    }
    ''')

    # Multiprocessing to use multiple processing cores.
    print(("Calculating thresholds (" + timeString() + "). This should only take a few minutes."))
    with ThreadPool(p) as pool:
        thresholds = pool.starmap(find_thresh, zip(itertools.repeat(get_thresh), matrices, itertools.repeat(pv),
                                  itertools.repeat(r_background)))

    print(("Total motifs read: " + str(len(thresholds))))
    print("Writing output file.")
    get_put_motifs(motif_file, motif_outfile, ow, thresholds)
    print(("Done (" + timeString() + ")"))

    return


if __name__ == '__main__':
    # Check if TFMPvalues package is installed, attempt to install if necessary.
    check_r_install('TFMPvalue')
    TFMPvalues = importr('TFMPvalue')
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-m", "--motif", dest="motif_file", required=True)
    parser.add_argument("-o", "--outfile", dest="motif_outfile", required=True)
    parser.add_argument("-a", "--a_freq", dest="a_freq", required=False, default=0.25, type=float)
    parser.add_argument("-c", "--c_freq", dest="c_freq", required=False, default=0.25, type=float)
    parser.add_argument("-g", "--g_freq", dest="g_freq", required=False, default=0.25, type=float)
    parser.add_argument("-t", "--t_freq", dest="t_freq", required=False, default=0.25, type=float)
    parser.add_argument("-p", "--processors", dest="processors", required=False, default=1, type=int)
    parser.add_argument("-pc", "--pseudocounts", dest="pseudocounts",
                        required=False, default=0.1, type=float)
    parser.add_argument("-pv", "--pval", dest="p_val",
                        required=False, default=0.00001, type=float)
    parser.add_argument("-ow", "--overwrite", action="store_true", required=False)
    args = parser.parse_args()

    if sum([args.a_freq, args.c_freq, args.g_freq, args.t_freq]) == 1:
        bp = [args.a_freq, args.c_freq, args.g_freq, args.t_freq]
    else:
        print("Background frequencies must equal 1. Check input parameters, exiting.")
        sys.exit()

    main(args.motif_file, args.motif_outfile, args.pseudocounts, bp, args.overwrite, args.p_val, args.processors)
