#!/usr/bin/env python3
"""
For the specified vcf and expression file, determine which genes are located within a specific range
from each variant. These are noted and their z-scores determined for samples with vs without the
variant.

Usage: python3 gene_expression.py -i <input.vcf> -e <expression.bed> -ov <output.vcf> -ob <output.bed> [OPTIONS]

Args:
    -i (str) = Name of vcf file to process.
    -e (str) = Name of expression file.
    -ov (str) = Name of vcf output file.
    -ob (str) = Name of bed output file that will be gene-based.
    -size (int, optional) = An integer to define the distance from the center or edges of the
        loci to look for gene overlap. 50 kb by default.
    -intergenic (bool, optional) = If included, will not print variants that directly overlap with a gene in
        the expression file to output. Good for regulatory element searches. False by default.
    -th (int, optional) = If given, will exclude genes that are below this threshold for all
        samples in the expression file. Good for excluding genes that aren't expressed. 0 by default.
"""
# TODO - This is very incomplete.
import argparse

# TODO - Scratch this class and use those from activity.py instead.
class Vcf:

    def __init__(self, sample_name, input_path, output_path):
        self.name = sample_name
        self.vcf_f = open(input_path)
        self.matches = []
        if output_path is not None:
            self.out_f = open(output_path, "w")
        else:
            self.out_f = None

        # Initialize current variant (aka current line).
        line = self.vcf_f.readline()
        info_needed = True
        info = ('##INFO=<ID=GZEXP,Number=.,Type=Float,Description="Gene z-score for Expression ',
                '(variant expression vs reference expression)">')

        # Skip info lines.
        while line.startswith("#"):
            # Print new info lines at the top of the ##INFO section.
            if self.out_f is not None:
                if info_needed and line.startswith("##INFO"):
                    print(info, file=self.out_f)
                    info_needed = False
                print(line, file=self.out_f, end="")
            line = self.vcf_f.readline()

        # Current variant.
        self.parse_line(line)

    # Outputs old variant if output path is given.
    def next_variant(self, options):
        if self.out_f is not None:
            self.output_var(options)

        self.matches = []
        line = self.vcf_f.readline()

        # skip info lines-
        while line.startswith("#"):
            line = self.vcf_f.readline()

        self.parse_line(line)

    # Remember to check if self.out_f is present before calling output_var.
    def output_var(self, options):

        # If the filter vcf option is on, don't print variants that don't
        # have an expression z-score (above the threshold if present).
        if options.filter_vcf and len(self.matches) == 0:
            return

        line_list = self.original_line.split("\t")

        # Initialize the line to be printed
        line = line_list[0]
        for idx in range(1, len(line_list)):

            if idx == 7 and len(self.matches) > 0:
                info = ""
                # Matches should have gene symbol, Position object, and zscore
                for match in self.matches:
                    (iden, pos, zscore) = match
                    if info != "":
                        info += ";"
                    info += "GZEXP=" + str(round(zscore, 3))
                if line_list[idx] != '.':
                    info += ';' + line_list[idx]
                line += '\t' + info
            else:
                line += '\t' + line_list[idx]

        print(line, file=self.out_f)

    def parse_line(self, line):
        """
        Read in a line of the vcf and create a Position object from the variant position.

        Args:
            line (str): VCF line to parse.

        Returns:
            None or a tuple with the following information (in order):
            The position of the variant as a Position class object (see below)
            A list of the motif names that matched that variant
            A list of other motif fields for futher processing
            A list of other INFO fields
            The original line
        """
        line = line.strip()

        line_list = line.split("\t")

        # Find variant position
        chrom = line_list[0]
        start_pos = int(line_list[1])
        ref_seq = line_list[3]
        end_pos = start_pos + len(ref_seq)
        var_pos = Position(chrom, start_pos, end_pos)

        self.var_pos = var_pos
        self.original_line = line

        return

    def close(self):
        self.vcf_f.close()

    def __str__(self):
        return self.name


class Position:
    """
    Args:
        chrom (string): Chromosome (chr1, chr2, etc).
        start = start position
        end = end position
    """

    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start = start_pos
        self.end = end_pos

    def overlaps(self, pos_b):
        """
        Determine whether this Position overlaps Position pos_b.

        Args:
            pos_b (Position): A Position object.

        Returns:
            (bool): Whether self overlaps with Position pos_b
        """
        if pos_b is None or self.chrom != pos_b.chr:
            return False

        start_max = max(self.start, pos_b.start)
        end_min = min(self.end, pos_b.end)

        return start_max <= end_min

    def add_wings(self, wing_length):
        """
        Add wings to a Position start and end and return the two
        resulting positions as a tuple.

        Args:
            position (int): Position to add the wings to.
            wing_length (int): The length of the wings to add to each side of the position.

        Returns:
            wing_positions (tuple): A tuple containing the position of each wing.
        """
        wing_length = int(wing_length)

        wing_start = self.start - wing_length
        wing_stop = self.end + wing_length

        wing_positions = (wing_start, wing_stop)

        return wing_positions

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)


def get_gene_exp(gene_f, thresh):
    """
    Grabs gene names and expression values for each sample.

    Args:
        gene_f (str): Expression file in bed-like format.
        thresh (float): Threshold used to filter genes that aren't expressed above this number in at least one sample.

    Returns:
        gene_dict (dict): Dictionary in following format
            {gene_name:(Position,[exp_vals])}
    """
    gene_dict = {}

    with open(gene_f) as f:

        for line in f:
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            pos = Position(chrom, start, stop)
            gene_symb = line[3]
            exp_vals = [float(x) for x in line[4:]]

            # Check if any of the data values are above the threshold. Skip if not.
            for data in exp_vals:
                if data >= thresh:
                    gene_dict[gene_symb] = (pos, exp_vals)
                    break

    return gene_dict


def main(inp_file, gene_file, wing_size, output, no_overlap=False, thresh=0):
    """
    Find expression differences for genes located within <wing_size> of each variant.

    Generate a z-score for expression of samples with the variant vs. without and report to output.

    Args:
        inp_file (str): Input file with locus positions and ID.
        gene_file (str): Gene annotation file in gtf format.
        wing_size (int): The length of the wings to add to each side the variant.
        output (str): Name of the output file.
        no_overlap (bool): Determines if loci that overlap the gene are included in output or not.
            Default=False (included).
        thresh (float): Genes that are below this threshold for all samples in the expression file will be excluded.
            Good for excluding genes that aren't expressed.
    """

    print("Parsing gene expression file.")
    gene_dict = get_gene_exp(gene_file, thresh)

    output_file = open(output, "w")

    with open(inp_file) as f:

        print("Comparing loci positions to gene positions.",
              " This will take a while for large datasets.")

        # Used to count the number with no genes that overlap.
        no_ovlp_count = 0
        ovlp_count = 0

        genes_ovlp = 0
        loci_excl_count = 0

        for line in f:

            # Create list to hold genes for that locus ID
            gene_list = []

            line = line.strip().split()

            # Assign line elements to variables
            locus_chrom = line[chrom_col]
            locus_start = int(line[start_col])
            locus_stop = int(line[end_col])

            # Check if the midpoint or edges of the locus are to be used
            if midpoint:
                # Get the midpoint of the locus
                locus_midpoint = Get_Midpoint(locus_start, locus_stop)

                # Add wings to the midpoint
                winged_positions = Add_Wings_Mid(locus_midpoint, wing_size)
            else:
                winged_positions = Add_Wings_Edge(locus_start, locus_stop, wing_size)

            wing_start, wing_stop = winged_positions

            # Iterate through each gene to see if it overlaps the winged locus positions
            for gene in gene_dictionary:
                gene_position = gene_dictionary[gene]
                gene_chrom = gene_position[0]

                # Check if it's even on the same chromosome as the locus, if not, skip to
                # the next gene
                if gene_chrom == locus_chrom:
                    gene_start_stop = gene_position[1]
                    gene_start = gene_start_stop[0]
                    gene_stop = gene_start_stop[1]

                    # Check that the original start/stop of the locus isn't wider than with
                    # the wings added, and if so, just use the original start/stop
                    if int(locus_start) < int(wing_start) and int(locus_stop) > int(wing_stop):
                        wing_start = locus_start
                        wing_stop = locus_stop

                    # Check for overlap of the two ranges and add the gene to the gene_list if
                    # this returns True.
                    if Overlap(wing_start, wing_stop, gene_start, gene_stop):
                        gene_list.append(gene)

            # Make sure all TSSs associated with the gene don't overlap the locus if
            # excl_promoter==True.
            if excl_proms:

                ovlp = False

                # Iterate through the prom_dictionary to check each TSS associated with the gene
                for item in gene_list:
                    for prom in prom_dictionary:
                        if item in prom:
                            prom_position = prom_dictionary[prom]
                            prom_start_stop = prom_position[1]
                            prom_start = prom_start_stop[0]
                            prom_stop = prom_start_stop[1]

                            # Check for overlap between the promoter and the locus, and if it
                            # occurs, just move to the next line
                            ovlp = Overlap(locus_start, locus_stop, prom_start, prom_stop)

                            if ovlp:

                                break

                    if ovlp:
                        break

                if not ovlp:
                    if len(gene_list) > 0:
                        ovlp_count += 1
                        genes_ovlp += len(gene_list)
                        line.append(";".join(gene_list))
                        print(*line[0:], sep="\t", file=output_file)
                    else:
                        no_ovlp_count += 1
                        line.append("NA")
                        print(*line[0:], sep="\t", file=output_file)
                else:
                    loci_excl_count += 1

            elif no_ovlp:

                overlap = False

                # Iterate through the gene_list to check for overlap
                for item in gene_list:

                    gene_position = gene_dictionary[item]
                    gene_chrom = gene_position[0]

                    # Check if it's even on the same chromosome as the locus, if not, skip to
                    # the next gene
                    if gene_chrom == locus_chrom:
                        gene_start_stop = gene_position[1]
                        gene_start = gene_start_stop[0]
                        gene_stop = gene_start_stop[1]

                        # Check for overlap between the promoter and the locus, and if it occurs,
                        # just move to the next line
                        overlap = overlap(locus_start, locus_stop, gene_start, gene_stop)

                    if overlap:
                        break

                if not overlap:
                    if len(gene_list) > 0:
                        ovlp_count += 1
                        genes_ovlp += len(gene_list)
                        line.append(";".join(gene_list))
                        print(*line[0:], sep="\t", file=output_file)
                    else:
                        no_ovlp_count += 1
                        line.append("NA")
                        print(*line[0:], sep="\t", file=output_file)
                else:
                    loci_excl_count += 1

            else:

                # Print to the output file
                # Check if there are genes in the gene_list, if not it will print NA
                if len(gene_list) > 0:
                    ovlp_count += 1
                    genes_ovlp += len(gene_list)
                    line.append(";".join(gene_list))
                    print(*line[0:], sep="\t", file=output_file)
                else:
                    no_ovlp_count += 1
                    line.append("NA")
                    print(*line[0:], sep="\t", file=output_file)

    if excl_proms:
        print("Number of loci excluded for overlap with a TSS: " + str(loci_excl_count))
    if no_overlap:
        print("Number of loci excluded for overlap with a gene: " + str(loci_excl_count))

    print("Number of loci with no overlapping genes: " + str(no_ovlp_count))
    print("Number of loci with overlapping genes: " + str(ovlp_count))

    output_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-g", "--gene", dest="gene_file", required=True)
    parser.add_argument("-size", "--wingsize", dest="wing_size", type=int, default=50000)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-noovlp", "--nooverlap", action="store_true")
    parser.add_argument("-th", "--threshold", dest="thresh", type=float, required=False, default=0)

    args = parser.parse_args()

    print("Input file {}\n Gene file {}\n Wing size {}\n Output {}\n\n".format(
        args.input_file,
        args.gene_file,
        args.wing_size,
        args.output_file
    ))

    # Easier to use argument variables
    inp_file = args.input_file
    g_file = args.gene_file
    wing_size = args.wing_size
    out_file = args.output_file
    no_ovlp = args.nooverlap
    thresh = args.thresh

    main(inp_file, g_file, wing_size, out_file, no_ovlp, thresh)
