#!/usr/bin/env python
"""
For a given motif annotated vcf file (already run through motifs.py),
remove all motif matches for TFs that are not expressed in at least one sample above the threshold.
    
Usage: tf_expression.py -i <input.vcf> -e <expression.bed> -o <output.txt> [OPTIONS]

Args:
    -i (required) <input.vcf>: Name of sorted variant file to process. 
    -o (required) <output.vcf>: Name of output file to be created.
    -e (required) <expression.bed>: An expression 'bed' file.
    -th (optional) <5>: TFs are considered expressed if they are above this threshold. 
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


def filter_motifs(motifs, gene_dict, thresh):
    """
    Get lists of motifs above threshold and below for a list of motif names & the gene dict.

    Args:
        motifs (list): Motif names for a given variant.
        gene_dict (dict): Dictionary containing gene names and expression values.
        thresh (float): Expression threshold that TFs must meet to be included in output.

    Returns:
        passed (list): List of motif indices that pass threshold.
        failed (list): List of motif indices that fail threshold.
    """

    passed = []
    failed = []

    for item in motifs:
        pass_th = False
        
        try:
            # TODO - Try to handle complexes and such - 'ATX::TCF3', etc.
            exp_vals = gene_dict[item]
        except:
            print(item + " not found in gene dict.")
            continue

        # Check if any of the expression values meet the threshold.
        for x in exp_vals:
            if x >= thresh:
                pass_th = True
                motif_idx = motifs.index(item)
                passed.append(motif_idx)
                break

        if pass_th == False:
            motif_idx = motifs.index(item)
            failed.append(motif_idx)

    return (passed, failed)


def process_line(line, output_f, gene_dict, thresh):
    """
    Process a vcf record and print to output.

    Args: 
        line (str): Line to parse.
        output_f (file object): Output file. 
        gene_dict (dict): Dict with gene names as keys and expression values for samples in vcf as
            values.
        thresh (float): Expression threshold that TFs must meet to be included in output.
    """
    line = line.strip()
    line_list = line.split("\t")
    info_fields = line_list[7].split(";")
    (motifns, motif_other, other_info) = ([], {}, [])

    for field in info_fields:
        if field.startswith("MOTIFN="):
            motifns = field[7:].split(',')
            # TODO - Check how lines with no motif matches are handled. Is MOTIFN still present? 
            # If not, skip line.

            # Get motifs that pass expression threshold.
            passed_idx, failed_idx = filter_motifs(motifns, gene_dict, thresh)  
            for i in sorted(failed_idx, reverse=True):
                del motifns[i]
        # Create dict from other motif fields.
        # TODO - Check if there is ALWAYS an entry for the 'motif_other' fields
        # corresponding to the indices in the 'failed_idx' list.
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

    new_info_fields = "MOTIFN=" + ",".join(motifns)

    # Delete matches for the other motif info fields that didn't pass the threshold.
    for f in motif_other:
        motif_other[f] = vals
        for i in sorted(failed_idx, reverse=True):
            del vals[i]
        motif_other[f] = vals
        new_info_fields += ";" + f + ",".join(vals)

    for x in other_info:
        new_info_fields += ";" + x

    line_list[7] = new_info_fields
    new_line = "\t".join(line_list)  # New, filtered output line.
    print(new_line, file=output_f)

    return 


def get_genes(exp_file, samples, threshold):
    """
    Reads in and parses the .bed expression file.

    Args: 
        exp_file (str): Name of expression file.
        samples (list): Names of the samples in the vcf file.
        threshold (float): Expression threshold to filter lowly/unexpressed genes.

    Returns: 
        gene_dict (dict): {gene_name: [expression_vals]}. Only include values for samples in 
            the vcf.
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


def main(inp_file, exp_file, out_file, th=5):
    """
        For a given motif annotated vcf file (already run through motifs.py), remove all motif 
        matches for TFs that are not expressed in at least one sample above the threshold.

        Args:
            -i (str): Name of sorted variant file to process. 
            -o (str): Name of output file to be created.
            -e (str): Name of expression file.
            -th (float): TFs are considered expressed if they are above this threshold. 
    """

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

    output_f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-e", "--expression", dest="exp_file", required=True)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold",
                        required=False, default=5)

    args = parser.parse_args()

    # Easier to use argument variables
    inp_file = args.input_file
    exp_file = args.exp_file
    out_file = args.output_file
    th = float(args.threshold)

    main(inp_file, exp_file, out_file, th)