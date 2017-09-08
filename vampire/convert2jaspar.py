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
        Input motif format. Options: 'transfac', 'encode', 'homer', 'minmeme'
        'MINMEME' are the .meme database files.
"""

import Bio.motifs as bmotifs
from Bio.motifs import jaspar
import argparse
import numpy as np
import sys


# def convert_meme(input_file, output_file):
#     """
#     Convert file of MEME motifs to JASPAR (2016) format.
#     """
#     with open(input_file) as f:
#         motifs = []
#         for m in bmotifs.parse(f, "MEME"):
#             motifs.append(bmotifs.write(m, "JASPAR"))

#     out_file = open(output_file, "w")
#     print("\n".join(motifs), file=out_file)
#     out_file.close()

#     return


def convert_min_meme(input_file, output_file):
    """
    Convert file of minimal MEME motifs to JASPAR (2016) format.
    """
    with open(input_file) as f:
        info = ""
        pfm = []
        first = True  # Keep track of first line.
        out_file = open(output_file, "w")

        for line in f:
            if line.startswith("MOTIF"):
                if not first:
                    print(info, file=out_file)
                    pfm = np.array(pfm)
                    out_pfm = np.transpose(pfm)  # Permute columns to rows.
                    out = ["\t".join(x) for x in out_pfm]
                    print("\n".join(out), file=out_file)
                    pfm = []  # Make new empty list.
                name_id = line.strip().split()
                name = name_id[2].split("_")[0]
                info = ">" + name_id[1] + " " + name.strip(")(")
                first = False
            # Skip empty and extraneous info lines.
            elif line.strip() and not line.startswith("letter") and line.strip()[0].isdigit():
                line = line.strip().split()
                pfm.append(line)

        print(info, file=out_file)
        pfm = np.array(pfm)
        out_pfm = np.transpose(pfm)  # Permute columns to rows.
        out = ["\t".join(x) for x in out_pfm]
        print("\n".join(out), file=out_file)

        out_file.close()

    return


def convert_transfac(input_file, output_file):
    """
    Convert file of TRANSFAC motifs to JASPAR (2016) format.
    """
    out = ""
    with open(input_file) as f:
        mtifs = bmotifs.parse(f, "TRANSFAC")
    for x in mtifs:
        try:
            x.name = x['NA']
        except KeyError:
            x.name = None

        try:
            x.matrix_id = x['ID']
        except KeyError:
            x.matrix_id = None

    out = jaspar.write(mtifs, 'jaspar')

    out_file = open(output_file, "w")
    print(out, file=out_file)
    out_file.close()

    return


def convert_encode(input_file, output_file):
    """
    Convert file of ENCODE pfm motifs to JASPAR (2016) format.
    """
    with open(input_file) as f:
        info = ""
        pfm = []
        first = True  # Keep track of first line.
        out_file = open(output_file, "w")

        for line in f:
            if line.startswith(">"):
                if not first:
                    print(info, file=out_file)
                    pfm = np.array(pfm)
                    out_pfm = np.transpose(pfm)  # Permute columns to rows.
                    out = ["\t".join(x) for x in out_pfm]
                    print("\n".join(out), file=out_file)
                    pfm = []  # Make new empty list.
                name_id = line.strip(">").split()
                info = ">" + name_id[1] + " " + name_id[0].split("_")[0]
                first = False
            elif line.strip():  # Skip empty lines.
                line = line.strip().split()[1:]  # Skip first character in line for dominant base at position.
                pfm.append(line)

        print(info, file=out_file)
        pfm = np.array(pfm)
        out_pfm = np.transpose(pfm)  # Permute columns to rows.
        out = ["\t".join(x) for x in out_pfm]
        print("\n".join(out), file=out_file)

        out_file.close()

    return


def convert_homer(input_file, output_file):
    """
    Convert file of HOMER pfm motifs to JASPAR (2016) format.
    """
    with open(input_file) as f:
        info = ""
        pfm = []
        first = True  # Keep track of first line.
        out_file = open(output_file, "w")

        for line in f:
            if line.startswith(">"):
                if not first:
                    print(info, file=out_file)
                    pfm = np.array(pfm)
                    out_pfm = np.transpose(pfm)  # Permute columns to rows.
                    out = ["\t".join(x) for x in out_pfm]
                    print("\n".join(out), file=out_file)
                    pfm = []  # Make new empty list.
                name_id = line.strip(">").split()
                info = ">" + name_id[1] + " " + name_id[0].split("_")[0]
                first = False
            elif line.strip():  # Skip empty lines.
                line = line.strip().split()[1:]  # Skip first character in line for dominant base at position.
                pfm.append(line)

        print(info, file=out_file)
        pfm = np.array(pfm)
        out_pfm = np.transpose(pfm)  # Permute columns to rows.
        out = ["\t".join(x) for x in out_pfm]
        print("\n".join(out), file=out_file)

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
            Input motif formats. Valid options: 'transfac', 'encode', 'homer', 'minmeme'
    """

    # if m_type == "MEME":
    #     print("Converting from MEME to JASPAR format.")
    #     convert_meme(input_file, output_file)
    if m_type == "MINMEME":
        print("Converting from minimal MEME to JASPAR format.")
        convert_min_meme(input_file, output_file)
    elif m_type == "TRANSFAC":
        print("Converting from TRANSFAC to JASPAR format.")
        convert_transfac(input_file, output_file)
    elif m_type == "ENCODE":
        print("Converting from ENCODE to JASPAR format.")
        convert_encode(input_file, output_file)
    elif m_type == "HOMER":
        print("Converting from HOMER to JASPAR format.")
        convert_homer(input_file, output_file)
    else:  # Should never be the case.
        print("Motif format not found, exiting.")
        sys.exit()

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_file", required=True)
    parser.add_argument("-m", "--motif", dest="m_type", action='store',
                        choices=['TRANSFAC', 'ENCODE', 'HOMER', 'MINMEME'], required=True, type=str.upper)

    args = parser.parse_args()

    main(args.inp_file, args.out_file, args.m_type)
