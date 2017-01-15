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

    # TODO - Implement bed output.
    -ob (str) = Name of bed output file that will be gene-based.

    -size (int, optional) = An integer to define the distance from the center or edges of the
        loci to look for gene overlap. 50 kb by default.

    # TODO - Implement this. Should really use an annotation file that has exons for each gene so that
    #   intronic REs are also included.
    -nc (bool, optional) = If included, will not print variants that directly overlap with a gene in
        the expression file to output. Good for regulatory element searches. False by default.

    -th (int, optional) = Z-score magnitude threshold that must be met for variants/loci to be reported to output.
        Default is 0, so all loci a variant overlaps will be reported. 0 by default.
    -eth (int, optional) = If given, will exclude genes that are below this threshold for all
        samples in the expression file. Good for excluding genes that aren't expressed. 0 by default.
"""
# TODO - This is very incomplete. # CCC-WK: Why?
from __future__ import print_function    # so Ninja IDE will stop complaining & show symbols
import argparse
import time
from statistics import mean, stdev


def timeString():
    """ Return time as a string YearMonthDay.Hour.Minute.Second
    """
    return time.strftime('%Y%m%d.%H:%M:%S')


# TODO - Move into a utils.py file and import as appropriate.
class Position(object):
    """
    Use to represent and handle genomic ranges more easily.

    Args:
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
    """

    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start = start_pos
        self.end = end_pos
        self.wing_start = None
        self.wing_end = None

    def overlaps(self, pos_b):
        """
        Return whether self overlaps Position pos_b's wings.

        Args:
            pos_b (Position): Another Position.

        Returns:
            bool: True if self overlaps with Position pos_b. False if not.
        """
        if pos_b is None:
            return False

        if self.chrom != pos_b.chrom:
            return False

        start1, start2, end1, end2 = (self.start, pos_b.wing_start, self.end, pos_b.wing_end)

        return end1 >= start2 and end2 >= start1

    def set_wings(self, wing_length):
        """
        Set wing positions for a Position object.

        Args:
            position (int): Position to add the wings to.
            wing_length (int): The length of the wings to add to each side of the position.

        Returns:
            wing_positions (tuple): A tuple containing the position of each wing.
        """
        wing_length = int(wing_length)

        self.wing_start = self.start - wing_length
        self.wing_end = self.end + wing_length

        return

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)


# TODO - Move into a utils.py file and import as appropriate. Add doc_string.
class Variant(object):
    """
    Use to process and handle variant records from a VCF more easily. Create from line of VCF file.
    """

    def __init__(self, line, all_sample_names):
        self.line_list = line.strip().split("\t")
        self.pos = Position(self.line_list[0], int(self.line_list[1]),
                            (int(self.line_list[1]) + len(self.line_list[3])))
        self.ref_allele = self.line_list[3]
        self.var_allele = self.line_list[4]
        self.iden = self.line_list[2]
        self.orig_line = line.strip()
        self.info_fields = self.line_list[7].split(";")
        self.var_samples, self.motif_fields = self.parse_info_fields()
        self.ref_samples = [x for x in all_sample_names if x not in self.var_samples]
        self.loci = []

        if self.var_samples is not None:  # Should never evaluate to False.
            self.num_var_samps = len(self.var_samples)
        else:
            self.num_var_samps = 0

    def parse_info_fields(self):
        """
        Get names of samples containing variant and motif INFO fields from a variant record's INFO fields.

        Args:
            self (Variant): Variant object.

        Returns:
            samples (list of str): List of samples in which variant was called.
            motif_fields (list of str): List of INFO fields for variant that contain MOTIF related information.
        """
        samples = None
        motif_fields = []

        for field in self.info_fields:
            if field != "INDEL":  # Take care of INDEL flag.
                field_info = field.split("=")

                # TODO - This is a hack work around a bug that's messing up the MOTIFN field in tf_expression.py.
                # Go back and actually figure out why the MOTIFN field is getting split up sometimes.
                try:
                    name, data = (field_info[0], field_info[1])
                except:
                    name, data = "BROKEN", None
            else:
                name, data = "INDEL", None

            # TODO - Write method that parses header to determine # samples with variant rather than this lazy method.
            if name == "SAMPSTV":
                samples = data.split(",")
            elif name.startswith("MOTIF"):
                motif_fields.append(field)

        return (samples, motif_fields)

    def get_variant_output(self, include_vcf=False):
        """
        Create VCF output line for given Variant object.

        Args:
            include_vcf (bool): True if variants that don't pass the z-score threshold for any Locus should excluded
                from output. False if they should be included.

        Returns:
            output (str): Line for Variant in appropriate VCF format.
                or
            None: If include_vcf is True and no Locus that Variant overlaps hits the z-score threshold.
        """
        info = self.info_fields
        info.insert(0, "EXPR=" + ",".join([x.ref_samples[self] for x in self.loci][0]))
        info.insert(0, "EXPV=" + ",".join([x.var_samples[self] for x in self.loci][0]))

        # TODO - Check and make sure next two lines are functioning properly.
        info.insert(0, "EXPNR=" + str(self.loci[0].num_valid_ref[self]))
        info.insert(0, "EXPNV=" + str(self.loci[0].num_valid_var[self]))

        # Use lists to maintain order in output so that GENE, EXPVZ, EXPTHN fields can all be matched up.
        z_scores = []
        pass_thresh = []
        loci_idens = []
        for item in self.loci:
            loci_idens.append(item.iden)
            pass_thresh.append(item.num_pass_thresh[self])
            tmp = "(" + ",".join([str(round(x, 4)) for x in item.z_scores[self][0]]) + ")"
            z_scores.append(tmp)
        info.insert(0, "EXPTHN=" + ",".join([str(x) for x in pass_thresh]))
        info.insert(0, "EXPVZ=" + ",".join(z_scores))
        info.insert(0, "GENE=" + ",".join(loci_idens))

        # Check if any loci have samples that pass the Z-score threshold.
        # TODO - Change this so it check that the NUMBER OF SAMPLES reaching the z-score threshold are enough.
        if any([x >= 1 for x in pass_thresh]):
            self.info_fields = info
            self.line_list[7] = ";".join(self.info_fields)
            output = "\t".join(self.line_list)
            return output
        else:
            return None


# TODO - Move into a utils.py file and import as appropriate. Make ActLocus a sub-class of Locus along with GeneLocus.
class Locus(object):
    """
    Use to process and handle loci records from an activity file more easily.

    Args:
        pos (Position): Position object holding genomic position of locus.
        orig_line (str): String from which the object was originally created.
        iden (str): Unique identifier for the locus.
        data (list of float): Data values for each sample for the record.
    """

    def __init__(self, line):
        line_list = line.strip().split("\t")
        self.pos = Position(line_list[0], int(line_list[1]), int(line_list[2]))
        self.orig_line = line.strip()
        self.iden = str(line_list[3])
        self.data = [float(x) for x in line_list[4:]]
        self.var_samples = {}
        self.ref_samples = {}
        self.ref_scores = {}
        self.var_scores = {}
        self.num_valid_ref = {}
        self.num_valid_var = {}
        self.num_pass_thresh = {}
        self.variants = []
        self.z_scores = {}

    def add_variant(self, variant, var_samples, ref_samples):
        """
        Add Variant object variant to list of Variants that overlap the Locus.
        """
        self.variants.append(variant)
        self.ref_scores[variant] = []
        self.var_scores[variant] = []
        self.var_samples[variant] = var_samples
        self.ref_samples[variant] = ref_samples
        self.num_valid_ref[variant] = len(ref_samples)
        self.num_valid_var[variant] = len(var_samples)
        self.num_pass_thresh[variant] = 0
        self.z_scores[variant] = []

    def calc_z_score(self, ref_ind, var_ind, variant, thresh=0):
        # TODO - Add docstring here.
        self.num_valid_ref[variant] = len(ref_ind)
        self.num_valid_var[variant] = len(var_ind)

        if len(ref_ind) > 1:  # If all samples (or all but 1) have the variant, can't calc z-score, return 'NA'.
            for entry in ref_ind:
                scores = self.ref_scores[variant]
                scores.append(self.data[int(entry)])
                self.ref_scores[variant] = scores
            for entry in var_ind:
                scores = self.var_scores[variant]
                scores.append(self.data[int(entry)])
                self.var_scores[variant] = scores

            ref_mean = mean(self.ref_scores[variant])
            ref_std = stdev(self.ref_scores[variant])

            if ref_std == 0:  # If only one sample has ref, will have no variance. Should never happen.
                for i in var_ind:
                    vals = self.z_scores[variant]
                    vals.append("NA")
                    self.z_scores[variant] = vals
            else:
                score = [((x - ref_mean) / ref_std) for x in self.var_scores[variant]]

                for item in score:
                    if abs(item) >= thresh:  # Check number of variant samples that passed the threshold.
                        passed = self.num_pass_thresh[variant]
                        passed += 1
                        self.num_pass_thresh[variant] = passed
                    vals = self.z_scores[variant]
                    vals.append(score)
                    self.z_scores[variant] = vals
        else:
            for i in var_ind:
                vals = self.z_scores[variant]
                vals.append("NA")
                self.z_scores[variant] = vals
        return


def get_expression_samples(header_line):
    """
    Parse header of expression file to return sample names and column indices.

    Args:
        header_line (str): Header line from expression file.

    Returns:
        act_samples (dict): Dictionary of {sample_name (str): sample_data_index (int)}.
            sample_data_index is index for data in sample list, not the line as a whole.
            e.g.: [samp1, samp2, samp3] & [20, 10, 5] for data values, then {'samp1': 0}.
    """
    line_list = header_line.strip().split("\t")
    samples = line_list[4:]
    exp_samples = {}

    for item in samples:
        samp_idx = samples.index(item)
        sample = item.split(".")[0]
        exp_samples[sample] = samp_idx

    return exp_samples


def get_vcf_samples(header_line):
    """
    Parse header of VCF file to return sample names found in file.

    Args:
        header_line (str): Header line from VCF file.

    Returns:
        vcf_samples (list of str): List of samples with data in VCF file.
    """
    line_list = header_line.strip().split("\t")
    samples = line_list[9:]
    vcf_samples = []

    for item in samples:
        sample = item.split(".")[-1]
        if sample not in vcf_samples:
            vcf_samples.append(sample)

    return vcf_samples


def compare_samples(exp_samples, vcf_samples):
    """
    Compare samples from expression file and vcf file.

    Return only samples in both as a list and delete those not found in both from the exp_samples dict.

    Args:
        exp_samples (dict): {(aexp_sample_names (str)): sample_indices (int)}
        vcf_samples (list of str): List of samples found in VCF file.

    Returns:
        common_samps (list of str): List of names of samples found in both the expression file and VCF file.
        valid_exp_samps (dict): Dict of {sample_names (str): expression file data column index (int)} for samples found
            in both the expression file and VCF file.
    """
    common_samps = list(set(list(exp_samples)) & set(vcf_samples))
    valid_exp_samples = {}

    # Create new dict for expression samples containing only those found in VCF file as well.
    for x in common_samps:
        valid_exp_samples[x] = exp_samples[x]

    return (common_samps, valid_exp_samples)


def parse_expression_file(expression_file, thresh):
    """
    Parse expression file to get data values for each record along with sample
    names and indices.

    Args:
        expression_file (str): Path to expression file to process.
        thresh (float): Genes with no values above this threshold will be excluded from analysis.

    Returns:
        exp_samples (dict): Dict of {sample_name: index for expression vals}.
        exp_data (list of Locus): List of Locus objects.
    """
    with open(expression_file) as f:
        header = f.readline().strip()
        exp_samples = get_expression_samples(header)  # Get sample names/indices.
        exp_data = []

        for line in f:
            vals = line.strip().split()[4:]
            if any(float(x) > thresh for x in vals):
                record = Locus(line)
                exp_data.append(record)

        return (exp_samples, exp_data)


def main(vcf_file, exp_file, out_vcf, thresh=0, size=50000, include_vcf=False, ethresh=0):
    """
    Compare expression of loci for samples harboring a variant within a given locus to those samples that do not.

    For a given motif annotated VCF file (already run through motifs.py) and a bed-like file for loci of interest and
    some value for each loci for each sample, find loci that overlap a variant and compare the value of samples with
    the variant to those without the variant. Report z-scores for each loci overlapped in an output VCF and report the
    variants for each loci in a bed-like, loci-centric output file as well.

    Args:
        vcf_file (str): Path to sorted variant file to process.
        exp_file (str): Path to expression 'bed' file.
        out_vcf (str): Path to VCF output file to be created.
        out_bed (str): Path to loci output file to be created.
        thresh (float, optional): Z-score magnitude that must be met for variants/loci to be reported to output.
        ethresh (float, optional): Expression threshold below which genes will be excluded from analysis.
        size (int, optional): Set number of samples that must meet z-score threshold for locus to be reported to
            bed output file. So this number of samples must have the variant and be significantly affected by it.
        include_vcf (bool, optional): True if variants should be reported in the VCF output even if they don't lie near
            a Locus and significantly affect its expression.
    """
    print("Parsing expression data file at: " + timeString() + ".")
    exp_samps, exp_data = parse_expression_file(exp_file, ethresh)

    # Set wing boundaries.
    for item in exp_data:
        item.pos.set_wings(size)

    output_vcf = open(out_vcf, "w")

    with open(vcf_file) as f:

        # Add new INFO lines.
        line = f.readline().strip()
        now = time.strftime("%c")
        info_needed = True
        # TODO - Refactor this so output isn't such an enormous mess. One info field, multiple sub-fields per motif.
        # TODO - Add sample names for those that pass threshold.
        info = '##INFO=<ID=GENE,Number=.,Type=String,Description="IDs for genes that variant overlaps.">\n'
        info += ('##INFO=<ID=EXPR,Number=.,Type=String,Description="Samples with the reference allele and '
                 'expression data.">\n')
        info += ('##INFO=<ID=EXPV,Number=.,Type=String,Description="Samples with the variant allele and expression '
                 'data.">\n')
        info += ('##INFO=<ID=EXPVZ,Number=.,Type=String,Description="Z-score for each gene near the variant.'
                 ' Calculated for each sample containing the variant for each gene.">\n')
        info += ('##INFO=<ID=EXPTHN,Number=.,Type=Integer,Description="Number of samples in which the variant meets'
                 'the z-score expression magnitude threshold.">\n')
        info += ('##INFO=<ID=EXPNV,Number=1,Type=Integer,Description="Number of samples containing variant'
                 ' and having expression data.">\n')
        info += ('##INFO=<ID=EXPNR,Number=1,Type=Integer,Description="Number of samples containing reference'
                 ' and having expression data.">')
        command = ('##venusaur=<ID=expression,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --expression ' + exp_file + ' --outputvcf ' + out_vcf +
                   ' --threshold ' + str(thresh) + ' --include_vcf ' + str(include_vcf) + '">')

        # Print new info lines at the top of the ##INFO section.
        while line.startswith("##"):
            if info_needed and line.startswith("##INFO"):
                print(command, file=output_vcf)
                print(info, file=output_vcf)
                info_needed = False
            print(line, file=output_vcf)
            line = f.readline().strip()

        vcf_samples = get_vcf_samples(line)  # Parse VCF sample header line to get samples present in file.
        print(line, file=output_vcf)

        print("Comparing samples in VCF file and expression file to find commonalities.\n")
        print("VCF samples: ", *vcf_samples, end="\n\n")
        print("Expression samples: ", *list(exp_samps.keys()), end="\n\n")

        common_samps, valid_exp_samps = compare_samples(exp_samps, vcf_samples)  # Get common samples b/twn the two.
        print("Common samples: ", *common_samps, end="\n\n")
        print("Processing variants. This may take some time.")
        # TODO - Progress bar might actually be a decent addition.

        for line in f:
            current_var = Variant(line, vcf_samples)
            loci_ovlp_var = []

            # Check if any of the variant samples actually have expression data as well, skip if not.
            for x in current_var.var_samples:
                if x in common_samps:
                    for item in exp_data:
                        if current_var.pos.chrom != item.pos.chrom:
                            continue
                        elif current_var.pos.overlaps(item.pos):
                            loci_ovlp_var.append(item)
                    break

            # If variant overlaps no loci, print to output only if include_vcf option used.
            if not loci_ovlp_var:
                if include_vcf:
                    print(line.strip(), file=output_vcf)
                    continue
                else:
                    continue

            # Get expression data indices for both samples with variant and without.
            var_exp_indices = [valid_exp_samps[x] for x in current_var.var_samples if x in valid_exp_samps]
            ref_exp_indices = [valid_exp_samps[x] for x in valid_exp_samps if x not in current_var.var_samples]

            # Calculate z-scores.
            for x, loc in enumerate(loci_ovlp_var):
                var_samples = [x for x in current_var.var_samples if x in valid_exp_samps]
                ref_samples = [x for x in valid_exp_samps if x not in current_var.var_samples]
                loc.add_variant(current_var, var_samples, ref_samples)  # Add Variant to Locus object.
                loc.calc_z_score(ref_exp_indices, var_exp_indices, current_var, thresh)
                current_var.loci.append(loc)  # Add Locus object to given Variant.
                loci_ovlp_var[x] = loc

            vcf_out_line = current_var.get_variant_output(include_vcf)

            if vcf_out_line is not None:
                print(vcf_out_line, file=output_vcf)
            elif include_vcf:
                print(line.strip(), file=output_vcf)

    print("Complete at: " + timeString() + ".")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-e", "--expression", dest="expression_file", required=True)
    parser.add_argument("-ov", "--outputvcf", dest="output_vcf", required=True)
    parser.add_argument("-s", "--size", dest="wing_size", required=False, default=50000, type=int)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=0, type=float)
    parser.add_argument("-eth", "--ethreshold", dest="ethreshold", required=False, default=0, type=float)
    parser.add_argument("-iv", "--include_vcf", action="store_true", required=False)

    args = parser.parse_args()

    inp_file = args.input_file
    exp_file = args.expression_file
    vcf_out = args.output_vcf
    size = args.wing_size
    th = args.threshold
    eth = args.ethreshold
    include_vcf = args.include_vcf

    main(inp_file, exp_file, vcf_out, th, size, include_vcf, eth)
