#!/usr/bin/env python
"""
For a given motif annotated vcf file (already run through motifs.py)
remove all motif matches for TFs that are not expressed in at least one sample.
    
Usage: tf_expression.py -i <input.vcf> -e <expression.bed> -o <output.txt> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process. 
    -o (required) <output.vcf>: Name of output file to be created.
    -e (required) <expression.bed>: An expression 'bed' file.
    -th (optional) <5>: TFs are considered expressed if they are above this threshold. 
    -fe (optional flag): If used, variants that do not match any motif for 
        an expressed protein will not be included in the output file.
"""
import sys
import argparse
from math import ceil
from math import log10


def parse_header(line):
# TODO - Get sample names and list of columns for each

def parse_line(line):
    """
    Parse a vcf record and return various info.

    Args: 
        opened_file (file object): An open .vcf file. 

    Returns: 
        motifsn (list): List of motif names.
        motif_other (list): List of lists containing other motif info fields.
            each other field has the INFO name as the first object
            ex: ["MOTIFV",["0.5", "0.24"]] -> MOTIFV=0.5;0.24
        other_info (list): List of other INFO fields.
        line (str): Original line.
    """
    line = line.strip()

    # Find list of motifs that match
    line_list = line.split("\t")
    info_fields = line_list[7].split(";")
    # Holds other fields (other than motif fields which may be modified)
    other_info = []
    (motifns, motif_other) = ([], [])

    for field in info_fields:
        if field.startswith("MOTIFN="):
            # Cut off 'MOTIFN=' then convert to an array
            motifns = field[7:].split(',')
        # Check for other motif fields such as:
        #(MOTIFN)   MOTIFV  MOTIFR  MOTIFC
        # MOTIFVH   MOTIFRH MOTIFVG MOTIFRG
        elif (field.startswith("MOTIFV=") or field.startswith("MOTIFR=") or
              field.startswith("MOTIFC=")):
            name = field[:7]
            values = field[7:].split(',')
            motif_other.append([name, values])
        elif (field.startswith("MOTIFVH=") or field.startswith("MOTIFRH=") or
              field.startswith("MOTIFVG=") or field.startswith("MOTIFRG=")):
            name = field[:8]
            values = field[8:].split(',')
            motif_other.append([name, values])
        else:
            other_info.append(field)

    return (motifns, motif_other, other_info, line)


def get_genes(opened_file, sample, threshold):
    """
    Reads in and parses the .bed expression file.

    Args: 
        opened_file = an already open input .bed file
        sample (str): Name of the sample in question.
        threshold (float): Expression threshold to filter lowly/unexpressed genes.

    Returns: 
        gene_dict = (str -> float) a dictionary of:
        gene names -> expression levels
    """

    # Get sample column.
    # TODO - Only get expression values for samples in VCF.

    headers = opened_file.readline().strip().split('\t')
    sample_idx = -1
    for idx in range(len(headers)):
        if headers[idx].upper() == sample.upper():
            sample_idx = idx
            break

    if sample_idx == -1:
        print("Err: Sample <" + sample + "> not found in input file header.")
        return []

    gene_names = []
    exp_values = []

    # Read in expression levels for each gene
    for line in opened_file:
        gene_name = line.strip().split('\t')[3].upper()
        # Only add the gene to the list if its expression is above the threshold.
        exp_level = float(line_list[sample_idx])
        if exp_level >= threshold:
            # Add gene names and expression values to respective lists
            gene_names.append(gene_name)
            exp_values.append(exp_level)

    return dict(zip(gene_names, exp_values))


def print_record(line_tup, output_f, filt):
    """
    Print a record to the output file.

    If filt is True, then the record won't be printed to output unless it 
    has a motif match to a TF that is expressed in one of the samples.

    Args:
        line_tup (tup): Info tuple with the following information
            (list) Motif names as strings.
            (str list list) list of lists of other motif info
            (str list) other info fields
            (str) Original line
        output_f (string): Name of output file.
        filt (bool): True if only records with a motif for a TF meeting the expression
            threshold should be output. False otherwise.
    """
    # Unpack info tuple.
    (motifns, motifos, motifes, oth_info, line) = line_tup

    columns = line.strip().split('\t')

    # Strings for motif info fields
    str_m_os = ["" for idx in range(len(motifos))]

    # Stored as:
    # ["MOTIFC", ['Y','N','Y','Y']] -> MOTIFC=Y,N,Y,Y
    names = ""
    explevels = ""

    for idx in range(len(motifns)):
        if idx != 0:
            names += ","
            for idy in range(len(str_m_os)):
                str_m_os[idy] += ","
            explevels += ","
        names += motifns[idx]
        for idy in range(len(str_m_os)):
            str_m_os[idy] += motifos[idy][1][idx]
        explevels += round(motifes[idx], 4)

    # If there are no matches, print the line unchanged or filter it out (return
    # without printing)
    if len(motifns) == 0 and filt == True:
        return

    outline = ""
    idx = 0

    for col in columns:
        if outline != "":
            outline += "\t"
        if idx == 7:
            if len(motifns) != 0:
                outline += "MOTIFN=" + names + ";MOTIFE=" + explevels
                for idy in range(len(motifos)):
                    outline += ";" + motifos[idy][0] + str_m_os[idy]
                for field in oth_info:
                    outline += ";" + field
            else:
                info = ""
                for field in oth_info:
                    if info != "":
                        info += ","
                    info += field
                outline += info
        else:
            outline += col
        idx += 1

    print(outline, file=output_f)

    return


def main(inp_file, exp_file, out_file, th, filt):
    # TODO - Get sample names.

    gene_dict = []
    with open(exp_file) as opened_ef:
        gene_dict = get_genes(opened_ef, sample_name, th)

    if len(gene_dict) == 0:
        print("Error, no genes above threshold found in expression file.")
        sys.exit()

    # Open output file.
    output_f = open(out_file, "w")

    with open(inp_file) as vcf:

        line = vcf.readline()

        # Skip info lines.
        while line.startswith("##"):
            print(line, file=output_f, end="")
            line = vcf.readline()

        # First non-## line is the header.
        # TODO - Parse header.
        print(line, file=output_f, end="")

        for line in vcf:
            current_var = parse_line(line)


            (motifns, motifos, oi, line) = current_var

            # Filter motifs by expression
            f_motifns = []

            # Copy names into filtered motifs list
            f_motifos = []
            for idx in range(len(motifos)):
                # Stored as:
                # ["MOTIFC", ['Y','N','Y','Y']] -> MOTIFC=Y,N,Y,Y
                info_name = motifos[idx][0]
                f_motifos.append([info_name, []])

            # Expression level of motifs
            f_motifes = []

            for idx in range(len(motifns)):
                motif_name = motifns[idx].upper()
                # Genes are only in the dictionary if their expression is above
                # the given threshold. Add them to the filtered lists if they are
                # expressed.

                if motif_name in gene_dict:
                    # Add expression level
                    f_motifes.append(gene_dict[motif_name])
                    # Add name and other data
                    f_motifns.append(motifns[idx])
                    for idy in range(len(motifos)):
                        f_motifos[idy][1].append(motifos[idy][1][idx])

            # Output filtered motifs
            tup = (f_motifns, f_motifos, f_motifes, oi, line)
            print_record(tup, output_f, options)

            current_var = get_next_var(vcf)

    output_f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-e", "--expression", dest="exp_file", required=True)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold",
                        required=False, default=5)
    parser.add_argument("-fe", "--filter", action="store_true", required=False)

    args = parser.parse_args()

    # Easier to use argument variables
    inp_file = args.input_file
    exp_file = args.exp_file
    out_file = args.output_file
    th = float(args.threshold)
    filt = arg.filter

    main(inp_file, exp_file, out_file, th, filt)