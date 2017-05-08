#!/usr/bin/env python3
"""
Convert various motif file formats to JASPAR format.

Usage: convert2jaspar.py -i <motifs.txt> -o <output.txt> [OPTIONS]

Args:
    -i (str):
        Path to sorted variant file to process.
    -o (str):
        Name of output file.
    -t (str):
        Input motif format. Options: 'meme', 'transfac', 'encode', 'homer'
"""
# XXX - Incomplete.

import Bio.motifs as bmotifs


def main():
    """

    """

    return

if __name__ == '__main__':
    # Check if TFMPvalues package is installed, attempt to install if necessary.
    utils.check_r_install('TFMPvalue')
    TFMPvalues = importr('TFMPvalue')

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_pre", required=True)
    parser.add_argument("-d", "--dthresh", dest="d_thresh", required=False, default=0, type=float)
    parser.add_argument("-m", "--motif", dest="motif_file", required=False, default=None)

    args = parser.parse_args()

    main(args.inp_file, args.out_pre, args.d_thresh, args.motif_file, bp, args.pseudocounts)