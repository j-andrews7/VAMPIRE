#!/usr/bin/env python
"""
Calculate motif thresholds for a given motif file and print to a new file.

Usage: thresholds.py -m <motifs.txt> -o <output.txt> [OPTIONS]

Args:
    -m (str): Input filename of a file containing PWMs.
    -o (str): Output filename.
    -a (float, optional) <0.25>: Background probability for A nucleotides. If not
        provided then all are assumed to be equally likely (all are 0.25).
    -t (float, optional) <0.25>: Background probability for T nucleotides.
    -c (float, optional) <0.25>: Background probability for C nucleotides.
    -g (float, optional) <0.25>: Background probability for G nucleotides.
    -p (int, optinal) <1>: Processor cores to utilize. Will decrease computation time linearly.
    -pc (float, optional) <0.1>: Pseudocounts value to be added to all positions of the motif frequency matrix
        before calculating the probability matrix.
    -pv (float, optional) <0.00001>: P-value to be used for defining thresholds for each motifs.
        I don't recommend changing this.
    -ow (flag, optional): OverWrite: If present, thresholds already present in the input file will be
        replaced in the output file.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import argparse
from multiprocessing.dummy import Pool as ThreadPool
import itertools
from pytfmpval import tfmp
from timeit import default_timer as timer

from Bio import motifs

from motif import get_put_motifs
from utils import timeString
# from memory_profiler import profile
# import pdb


# @profile
def find_thresh(combo):
    """
    Calculate the detection thresholds for each pwm matrix given the nucleotide
    background frequencies and desired p-value.

    Args:
        combo (tuple): Three element tuple as defined below.
            matrix (tuple): The PWM for which a detection threshold will be found.
                (motif name, rows of matrix concatenated to single string with spaces between positions).
            pval (float): P-value to which threshold will be calculated.
            bg (list of floats): List containing background nucleotide frequencies [A, C, G, T]
    """

    matrix, pval, bg = (combo[0], combo[1], combo[2])

    start = timer()
    mat = tfmp.read_matrix(matrix[1], bg=bg, mat_type="pwm")
    thresh = tfmp.pval2score(mat, pval)
    # pdb.set_trace()
    del mat
    end = timer()

    print(matrix[0] + ": " + str(end - start))

    return (matrix[0], thresh)


def main(motif_file, motif_outfile, pc, bp, ow, pv, p):
    matrices = []
    background = {'A': bp[0], 'C': bp[1], 'G': bp[2], 'T': bp[3]}
    print(("Baseline nucleotide frequencies:\n\t" + str(background)))

    # Calculate thresholds using pytfmpval.
    print(("Reading in motifs."))
    fh = open(motif_file)
    for m in motifs.parse(fh, "jaspar"):
        pfm = m.counts.normalize(pseudocounts=pc)    # Create frequency matrix.
        pwm = pfm.log_odds(background)              # Calculate to log likelihoods vs background.

        # Create matrix string from motif pwm.
        mat = pwm[0] + pwm[1] + pwm[2] + pwm[3]
        mat = [str(x) for x in mat]
        mat = " ".join(mat)
        matrices.append((m.name, mat))

    fh.close()

    # Multiprocessing to use multiple processing cores.
    print(("Calculating thresholds (" + timeString() + "). This should only take a few minutes."))
    thresholds = []

    with ThreadPool(p) as pool:
        for x in pool.imap_unordered(find_thresh, zip(matrices, itertools.repeat(pv),
                                                      itertools.repeat(bp)), chunksize=16):
            thresholds.append(x)

    print(("Total motifs read: " + str(len(thresholds))))
    print("Writing output file.")
    get_put_motifs(motif_file, motif_outfile, ow, dict(thresholds))
    print(("Done (" + timeString() + ")"))

    return


if __name__ == '__main__':
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

    if sum([args.a_freq, args.c_freq, args.g_freq, args.t_freq]) == 1.0:
        bp = [args.a_freq, args.c_freq, args.g_freq, args.t_freq]
    else:
        print("Background frequencies must equal 1. Check input parameters, exiting.")
        sys.exit()

    main(args.motif_file, args.motif_outfile, args.pseudocounts, bp, args.overwrite, args.p_val, args.processors)
