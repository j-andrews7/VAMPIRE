#!/usr/bin/env python3
"""
Convert various motif file formats to JASPAR format.

Usage: convert2jaspar.py -i <motifs.txt> -o <output.txt> -m <input motif format>

Args:
    -i (str):
        Path to motif file to process.
    -o (str):
        Name of output file.
    -t (str):
        Input motif format. Options: 'meme', 'transfac', 'encode', 'homer'
"""
# XXX - Incomplete.

import Bio.motifs as bmotifs
import argparse
import sys


def convert_meme(input_file, output_file):
    """
    Convert file of MEME motifs to JASPAR (2016) format.
    """
    print("Converting from MEME to JASPAR format.")
    with open(input_file) as f:
        motifs = []
        for m in bmotifs.parse(f, "MEME"):
            motifs.append(bmotifs.jaspar.write(m, "JASPAR"))

    out_file = open(output_file, "w")
    print("\n".join(motifs), file=out_file)
    out_file.close()

    return


def convert_transfac(input_file, output_file):
    """
    Convert file of TRANSFAC motifs to JASPAR (2016) format.
    """
    print("Converting from TRANSFAC to JASPAR format.")
    with open(input_file) as f:
        motifs = []
        for m in bmotifs.parse(f, "TRANSFAC"):
            motifs.append(bmotifs.jaspar.write(m, "JASPAR"))

    out_file = open(output_file, "w")
    print("\n".join(motifs), file=out_file)
    out_file.close()

    return
    


def main(input_file, output_file, m_type):
    """
    Convert motif file of given format to JASPAR (2016) format.

    Args:
        input_file (str):
            Name of input motif file.
        output_file (str):
            Name of output motif file (in JASPAR format).
        m_type (str):
            Input motif formats. Valid options: 'meme', 'transfac', 'encode', 'homer'
    """

    if m_type is "MEME":
        print("Converting from MEME to JASPAR format.")
        convert_meme(input_file, output_file)


    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_file", required=True)
    parser.add_argument("-m", "--motif", dest="m_type", action='store',
                        choices=['MEME', 'TRANSFAC', 'ENCODE', 'HOMER'], required=True, type=str.upper)

    args = parser.parse_args()

    main(args.inp_file, args.out_file, args.m_type)
