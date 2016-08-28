#!/usr/bin/env python
"""
For a given motif annotated vcf file (already run through motifs.py)
remove all motif matches for TFs that are not expressed in at least one sample above the threshold.
    
Usage: tf_expression.py -i <input.vcf> -e <expression.bed> -o <output.txt> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process. 
    -o (required) <output.vcf>: Name of output file to be created.
    -e (required) <expression.bed>: An expression 'bed' file.
    -th (optional) <5>: TFs are considered expressed if they are above this threshold. 
    -fe (optional) <True>: If True, variants that do not match any motif for an expressed protein 
        will not be included in the output file. Set to False to print original line for variants 
        that lie in motifs for TFs that are expressed above the threshold.
"""
import sys
import argparse


def parse_header(line):
    """
    Get all sample names present in the vcf header. 

    Assume each column is in format 'extra_info.samplename'.

    Args:
        line (str): Header line.

    Returns:
        samples (list): List of sample names.
    """
    line = line.strip().split("\t")[9:]
    samples = []

    for x in line:
        samp = x.split('.')[-1]
        if samp not in samples:
            samples.append(samp)

    return samples


# TODO - Function to get lists of motifs above threshold and below for a list of motif names & 
# the gene dict.

def process_line(line, filt, output_f, gene_dict):
    """
    Process a vcf record and print to output.

    Args: 
        line (str): Line to parse.
        filt (bool): True if records that don't have a motif for a TF that's expressed
            should be excluded from output. False if they should still be included.
        output_f (file object): Output file. 
        gene_dict (dict): Dict with gene names as keys and expression values for samples in vcf as
            values.
    """
    line = line.strip()
    line_list = line.split("\t")
    info_fields = line_list[7].split(";")
    (motifns, motif_other, other_info) = ([], {}, [])

    for field in info_fields:
        if field.startswith("MOTIFN="):
            # Cut off 'MOTIFN=' then create a list.
            motifns = field[7:].split(',')
        # Create dict from other motif fields.
        elif (field.startswith("MOTIFV=") or field.startswith("MOTIFR=") or
              field.startswith("MOTIFC=")):
            name = field[:7]
            values = field[7:].split(',')
            motif_other[name] = values
        elif (field.startswith("MOTIFVH=") or field.startswith("MOTIFRH=") or
              field.startswith("MOTIFVG=") or field.startswith("MOTIFRG=")):
            name = field[:8]
            values = field[8:].split(',')
            motif_other[name] = values
        else:
            other_info.append(field)

    return 


def get_genes(exp_file, samples, threshold):
    """
    Reads in and parses the .bed expression file.

    Args: 
        exp_file (str): Name of expression file.
        samples (list): Names of the samples in the vcf file.
        threshold (float): Expression threshold to filter lowly/unexpressed genes.

    Returns: 
        gene_dict (dict): {gene_name: [expression_vals]. Only include values for samples in the vcf.
    """
    data_cols = []  # Used to hold indices for columns for which to get expression data.
    gene_dict = {}

    with open exp_file as f:
        header = f.readline().strip().split('\t')
        for samp in header[4:]:
            if samp in samples:
                data_idx = header.index(samp)
                data_cols.append(data_idx)

        # Read in expression levels for each gene.
        for line in f:
            line = line.strip().split('\t')
            gene_name = line[3].upper()
            exp_vals = []
            for idx in data_cols:
                exp_vals.append(line[idx])

            passed_filter = [x for x >= threshold in exp_vals]
            if len(passed_filter) > 0:
                gene_dict[gene_name] = exp_vals

        return gene_dict


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
                outline += "MOTIFN=" + names
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


    output_f = open(out_file, "w")

    with open(inp_file) as vcf:

        line = vcf.readline().strip()

        # Skip info lines.
        while line.startswith("##"):
            print(line, file=output_f)
            line = vcf.readline().strip()

        # First non-## line is the header. Get sample names and print to output.
        samples = parse_header(line)
        print(line, file=output_f)

        gene_dict = get_genes(exp_file, sample_names, th)

        if len(gene_dict) == 0:
            print("Error, no genes above threshold found in expression file.",
                "\nTry lowering the threshold and ensure the expression file",
                " has values in the range that you expect.")
            sys.exit()

        for line in vcf:
            current_var = process_line(line, filt, output_f)

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
            tup = (f_motifns, f_motifos, oi, line)
            print_record(tup, output_f, filt)

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
