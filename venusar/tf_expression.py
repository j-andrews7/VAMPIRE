#!/usr/bin/env python
"""
Input type 1: (-i, -o): For a given motif annotated vcf file (already run through motifs.py)
Input type 2: (-i, -m, -mo): For a given vcf file (used for samples), filter the motif file

With motifs from either input type and a gene expression file (-e), process the motif set to
remove all motif matches for TFs that are not expressed in at least one sample above the
specified threshold.

Usage: tf_expression.py -i <input.vcf> -e <expression.bed> -o <output.vcf> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process.
    -o (required) <output.vcf>: Name of output file to be created.
        Not created if using -m and -mo
    -e (required) <expression.bed>: An expression 'bed' file.
    -m (optional) <motif.txt>: Tab-delimited key file containing a frequency
        matrix with each row corresponding to a base and each column
        corresponding to a position (JASPAR format).
        If specified ignores input.vcf and output.vcf
    -mo (optional) <motif_output.txt>: Name of output motif file to be created.
        if blank and -m then creates .tf_filtered version of motif.txt file
    -th (optional) <5>: TFs are considered expressed if they are above this threshold.
"""
import sys
import argparse
import time
import motif


def parse_header(line):
    """
    Get all sample names present in the vcf header.

    Assume each column is in format 'extra_info.samplename'.

    Args:
        line (str): Header line.

    Returns:
        samples (list of str): List of sample names.
    """
    line = line.strip().split("\t")[9:]
    samples = []

    for x in line:
        samp = x.split('.')[-1]
        if samp not in samples:
            samples.append(samp)

    return samples


def filter_motifs(motifs, gene_dict, thresh):
    """
    Get lists of motifs above threshold and below for a list of motif names & the gene dict.

    Args:
        motifs (list of str): Motif names for a given variant.
        gene_dict (dict): Dictionary containing gene names and expression values.
        thresh (float): Expression threshold that TFs must meet to be included in output.

    Returns:
        passed (list of int): List of motif indices that pass threshold.
        failed (list of int): List of motif indices that fail threshold.
    """

    passed = []
    failed = []

    for item in motifs:
        pass_th = False

        try:
            # TODO - Try to handle complexes and such - 'ATX::TCF3', etc.
            exp_vals = gene_dict[item]    # this is the line that must be in try

            # Check if any of the expression values meet the threshold.
            for x in exp_vals:
                if float(x) >= thresh:
                    pass_th = True
                    indices = [i for i, x in enumerate(motifs) if x == item]  # Find multiple motifs for same TF.
                    for x in indices:
                        if x not in passed:
                            passed.append(x)
                    break    # breaks out of x in exp_vals

        except:
            # TODO - Add option to retain motifs that don't match a gene in the expression file.
            pass_th = False

        if pass_th is False:
            indices = [i for i, x in enumerate(motifs) if x == item]  # Find multiple motifs for same TF.
            for x in indices:
                if x not in failed:
                    failed.append(x)

    return (passed, failed)


def process_line(line, gene_dict, thresh):
    """
    Process a vcf record and return modified version of input line for output
    Returning line rather than direct print enable later integration with
    process stack class elements.

    Args:
        line (str): Line to parse.
        gene_dict (dict): Dict with gene names as keys and expression values for samples in vcf as values.
        thresh (float): Expression threshold that TFs must meet to be included in output.
    Returns:
        modified version of input line
    """
    line = line.strip()
    line_list = line.split("\t")
    info_fields = line_list[7].split(";")
    (motifns, motif_other, other_info) = ([], {}, [])

    motifn_present = False  # Skip line if no motif was matched.

    for field in info_fields:
        if field.startswith("MOTIFN="):
            motifn_present = True
            motifns = field[7:].split(',')    # the set of names

            # Get motifs that pass expression threshold.
            passed_idx, failed_idx = filter_motifs(motifns, gene_dict, thresh)

            # Delete by high index first so that lower indices aren't changed.
            for i in sorted(failed_idx, reverse=True):
                del motifns[i]

        # Create dict from other motif fields.
        elif (field.startswith("MOTIFV=") or field.startswith("MOTIFR=") or field.startswith("MOTIFC=")):
            name = field[:7]
            values = field[7:].split(',')
            motif_other[name] = values
        elif (field.startswith("MOTIFVH=") or field.startswith("MOTIFRH=") or field.startswith("MOTIFVG=") or
              field.startswith("MOTIFRG=")):
            name = field[:8]
            values = field[8:].split(',')
            motif_other[name] = values
        else:
            other_info.append(field)

    # Just go to next line if motif info isn't present.
    if motifn_present is False:
        return

    new_info_fields = "MOTIFN=" + ",".join(motifns)

    # Delete matches for the motif info fields that didn't pass the threshold.
    for f in motif_other:
        vals = motif_other[f]
        for i in sorted(failed_idx, reverse=True):
            del vals[i]
        motif_other[f] = vals
        new_info_fields += ";" + f + ",".join(vals)

    for x in other_info:
        new_info_fields += ";" + x

    line_list[7] = new_info_fields
    new_line = "\t".join(line_list)  # New, filtered output line.

    return new_line


def get_genes(exp_file, samples, threshold, max_only):
    """
    Reads in and parses the .bed expression file.
    File format expected to be:
        Whose format is tab seperated columns with header line:
        CHR  START  STOP  GENE  <sample 1>  <sample 2>  ...  <sample n>

    Args:
        exp_file (str): Name of expression file.
        samples (list): Names of the samples in the vcf file.
        threshold (float): Expression threshold to filter lowly/unexpressed genes.
        max_only (bool): if true, gene_dict value is 1 value = max expression
            if false gene_dict value is list of expression values
            YYY: WARNING: if want list to have meaning
                then values needs to be tied to header sample names

    Returns:
        gene_dict (dict): {gene_name: [expression_vals]}.
            Only include values for samples in the vcf.
    """
    data_cols = []
    gene_dict = {}

    print('start read exp_file:' + format(exp_file))

    if max_only:
        # read and only return max exp value in gene_dict
        with open(exp_file) as f:
            header = f.readline().strip().split('\t')
            for samp in header[4:]:
                if samp in samples:
                    data_idx = header.index(samp)
                    data_cols.append(data_idx)

            # Read in expression levels for each gene.
            for line in f:
                line = line.strip().split('\t')
                gene_name = line[3].upper()
                exp_val = -1e1000
                for idx in data_cols:
                    if float(line[idx]) > exp_val:
                        exp_val = float(line[idx])

                gene_dict[gene_name] = exp_val

    else:
        # read and return exp value list in gene_dict
        with open(exp_file) as f:
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

                gene_dict[gene_name] = exp_vals

    return gene_dict


def main(inp_file, exp_file, out_file, th=5, motif_file="", motifout_file="", use_vcf=True):
    """
        If use_vcf true then:
            For a given motif annotated vcf file (already run through motifs.py),
            remove all motif matches for TFs that are
            not expressed in at least one sample above the threshold.
        Else:
            For a given vcf file read the samples,
            then filter the motif_file, only outputing to motifout_file those items
            not expressed in at least one sample above the threshold

        Args:
            -i (str): Name of sorted variant file to process.
            -o (str): Name of output file to be created.
            -e (str): Name of expression file.
            -th (float): TFs are considered expressed if they are above this threshold.
    """

    if use_vcf:
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

            print("Creating gene dictionary for expression data.")
            gene_dict = get_genes(exp_file, samples, th, True)

            if len(gene_dict) == 0:
                print("Error, no genes above threshold found in expression file.",
                      "\nTry lowering the threshold and ensure the expression fil",
                      "e has values in the range that you expect.")
                sys.exit()

            print("Filtering motif info for TFs that don't meet the expression threshold of " +
                str(th) + ". Found " + format(len(gene_dict)) + " genes. Start processing vcf file.")
            for line in vcf:
                new_line = process_line(line, gene_dict, th)
                if new_line is not None:
                    print(new_line, file=output_f)

        output_f.close()
    else:
        # this version processes the motif file

        with open(inp_file) as vcf:

            line = vcf.readline().strip()

            # Skip info lines.
            while line.startswith("##"):
                line = vcf.readline().strip()

            # First non-## line is the header. Get sample names and print to output.
            samples = parse_header(line)

        # done with vcf file; only used to get samples

        print("Creating gene dictionary for expression data.")
        gene_dict = get_genes(exp_file, samples, th, True)

        if len(gene_dict) == 0:
            print("Error, no genes above threshold found in expression file.",
                  "\nTry lowering the threshold and ensure the expression fil",
                  "e has values in the range that you expect.")
            sys.exit()

        print("Filtering motif info for TFs that don't meet the expression threshold of " +
            str(th) + ". Found " + format(len(gene_dict)) + " genes. Start filtering motifs.")

        motif.get_filterbygene_put_motifs(motif_file, motifout_file, th, gene_dict)

    print("COMPLETE.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-e", "--expression", dest="exp_file", required=True)
    parser.add_argument("-m", "--motif_input", dest="motif_file", required=False)
    parser.add_argument("-mo", "--motif_output", dest="motif_out_file", required=False)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=5, type=float)

    args = parser.parse_args()

    inp_file = args.input_file
    exp_file = args.exp_file
    out_file = args.output_file
    th = args.threshold

    if args.motif_file is not None:
        motif_file = args.motif_file
        if args.motif_out_file is not None:
            motifout_file = args.motif_out_file
        else:
            motifout_file = motif_file.replace('.txt', '.tf_filtered.txt')
        main(inp_file, exp_file, out_file, th, motif_file, motifout_file, False)
    else:
        main(inp_file, exp_file, out_file, th, None, None, True)
