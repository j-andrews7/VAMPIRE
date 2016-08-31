#!/usr/bin/env python
"""
For a given motif annotated VCF file (already run through motifs.py) and a bed-like file for loci of interest and some
value for each loci for each sample, find loci that overlap a variant and compare the value of samples with the variant
to those without the variant. Report z-scores for each loci overlapped in an output VCF and report the variants for
each loci in a bed-like, loci-centric output file as well.

Usage: activity.py -i <input.vcf> -a <activity.bed> -ov <output.vcf> -ob <output.bed> [OPTIONS]

Args:
    -i (str): Path to sorted variant file to process.
    -a (str): Path to activity 'bed' file.
    -ov (str): Path to VCF output file to be created.
    -ob (str): Path to loci output file to be created.
    -th (float, optional): Z-score magnitutde threshold that must be met for variants/loci to be reported to output.
    -thresh (float, optional): Z-score magnitude that must be met for variants/loci to be reported to output.
        filter_num (int, optional): Set number of samples that must meet z-score threshold for locus to be reported to
            bed output file. So this number of samples must have the variant and be significantly affected by it.
        include_bed (bool, optional):
"""
import argparse
from statistics import mean, stdev

# TODO - Option to retain variants that lie outside these loci in VCF output.
# TODO - Option to retain loci that don't contain a variant in loci output.


class Position:
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

    def is_before(self, pos_b):
        """
        Returns true if self comes before Position pos_b without overlapping.

        Args:
            pos_b (Position): Another Position.

        Returns:
            bool: True if self overlaps with Position pos_b. False if not.
        """
        if pos_b is None:
            return True

        if self.chr == pos_b.chr:
            return self.end < pos_b.start

        return self.chr < pos_b.chr

    def overlaps(self, pos_b):
        """
        Return whether self overlaps Position pos_b.

        Args:
            pos_b (Position): Another Position.

        Returns:
            bool: True if self overlaps with Position pos_b. False if not.
        """
        if pos_b is None:
            return False

        if self.chrom != pos_b.chr:
            return False

        start_max = max(self.start, pos_b.start)
        end_min = min(self.end, pos_b.end)

        return start_max <= end_min


# TODO - Move into utils.py file and import as appropriate.
class Variant:
    """
    Use to process and handle variant records from a VCF more easily.
    """

    def __init__(self, line):
        line_list = line.strip().split("\t")
        self.pos = Position(line_list[0], int(line_list[1]), len(line_list[3]))
        self.orig_line = line.strip()
        self.info_fields = line_list[7].split(";")
        self.samples, self.motif_fields = self.parse_info_fields()

        if self.samples is not None:
            self.num_samps = len(self.samples)

    def parse_info_fields(self):
        """
        Get sample names and motif INFO fields from a variant record's INFO fields.

        Args:
            self (Variant): Variant object.

        Returns:
            samples (list of str): List of samples in which variant was called.
            motif_fields (list of str): List of INFO fields for variant that contain MOTIF related information.
        """
        samples = None
        motif_fields = []
        for field in self.info_fields:
            field_info = field.split("=")
            name, data = (field_info[0], field_info[1])

            # TODO - Include method that parses header to determine # samples rather than this lazy method.
            if name is "set":
                samples = data.split("-")
            elif name.startswith("MOTIF"):
                motif_fields.append(field)

        return (samples, motif_fields)


# TODO - Move into utils.py file and import as appropriate.
class Locus:
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
        self.pos = Position(line_list[0], int(line_list[1]), len(line_list[2]))
        self.orig_line = line.strip()
        self.iden = str(line[3])
        self.data = [float(x) for x in line[4:]]

    def calc_z_score(self, ref_ind, var_ind):
        self.num_ref = len(ref_ind)
        self.num_var = len(var_ind)
        ref_scores = []
        var_scores = []

        if len(num_ref) > 1:  # If all samples (or all but 1) have the variant, can't calc z-score. Just return 'NA'.
            for entry in ref_ind:
                ref_scores.append(data[int(entry)])
            for entry in var_ind:
                vat_scores.append(data[int(entry)])

            ref_mean = mean(ref_scores)
            ref_std = stdev(ref_scores)

            if ref_std != 0:  # If only one sample has ref, will have no variance.




def get_activity_samples(header_line):
    """
    Parse header of activity file to return sample names and column indices.

    Args:
        header_line (str): Header line from activity file.

    Returns:
        act_samples (dict): Dictionary of {sample_name (str): data_col_index (int)}.
    """
    line_list = header_line.strip().split("\t")
    samples = line_list[4:]
    act_samples = {}

    for item in samples:
        samp_idx = line_list.index(item)
        sample = item.split(".")[0]
        act_samples[sample] = samp_idx

    return act_samples


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
        samp_idx = line_list.index(item)
        sample = item.split(".")[-1]
        vcf_samples.append(sample)

    return vcf_samples


def compare_samples(act_samples, vcf_samples):
    """
    Compare samples from activity file and vcf file.

    Return only samples in both as a list and delete those not found in both from the act_samples dict.

    Args:
        act_samples (dict): {(act_sample_names (str)): sample_indices (int)}
        vcf_samples (list of str): List of samples found in VCF file.

    Returns:
        common_samps (list of str): List of names of samples found in both the activity file and VCF file.
        valid_act_samps (dict): Dict of {sample_names (str): activity file data column index (int)} for samples found
            in both the activity file and VCF file.
    """
    common_samps = list(set(list(act_samples)) & set(vcf_samples))
    valid_act_samples = {}

    # Create new dict for activity samples containing only those found in VCF file as well.
    for x in common_samps:
        valid_act_samples[x] = act_samples[x]

    return (common_samps, valid_act_samples)


def parse_activity_file(activity_file):
    """
    Parse activity file to get data values for each record along with sample
    names and indices.

    Args:
        activity_file (str): Path to activity file to process.

    Returns:
        act_samples (dict): Dict of {sample_name: index for activity vals}.
        act_data (list of Locus): List of Locus objects.
    """
    with open(activity_file) as f:
        header = f.readline.strip()
        act_samples = get_activity_samples(header)  # Get sample names/indices.
        act_data = []

        for line in f:
            record = Locus(line)
            act_data.append(record)

        return (act_samples, act_data)


def process_variant(line, act_samples, act_data, act_vars):
    """
    Process a variant record and print to output.

    Determine if the variant lies within one of the loci from the activity
    file, and if so, calculate a z-score for samples with the variant vs.
    those without it using the activity data. Print the variant to the
    output VCF file with the loci IDs and z-scores as new INFO fields.

    Args:
        line (str): Variant record to process from the input VCF file.
        act_samples (dict): Dict of {sample_name: index for activity vals}.
        act_data (dict): Dict of {(record Position, record iden): activity vals for record}.
        act_vars (dict): Dict of {(record Position, record iden): (Variants)}
    """


def main(vcf_file, act_file, out_vcf, out_bed, thresh=0, filter_num=0, include_bed=False, include_vcf=False):
    """
    Compare activity of loci for samples harboring a variant without them and those that do not.

    For a given motif annotated VCF file (already run through motifs.py) and a bed-like file for loci of interest and
    some value for each loci for each sample, find loci that overlap a variant and compare the value of samples with
    the variant to those without the variant. Report z-scores for each loci overlapped in an output VCF and report the
    variants for each loci in a bed-like, loci-centric output file as well.

    Args:
        vcf_file (str): Path to sorted variant file to process.
        act_file (str): Path to activity 'bed' file.
        out_vcf (str): Path to VCF output file to be created.
        out_bed (str): Path to loci output file to be created.
        thresh (float, optional): Z-score magnitude that must be met for variants/loci to be reported to output.
        filter_num (int, optional): Set number of samples that must meet z-score threshold for locus to be reported to
            bed output file. So this number of samples must have the variant and be significantly affected by it.
        include_bed (bool, optional): True if loci should be reported in the bed output even if they don't have a
            variant in them that significantly affects their activity.
        include_vcf (bool, optional): True if variants should be reported in the VCF output even if they don't lie in
            a Locus and significantly affect its activity.
    """

    act_samps, act_data = parse_activity_file(act_file)

    output_vcf = open(out_vcf, "w")
    output_bed = open(out_bed, "w")

    with open(vcf_file) as f:

        line = vcf.readline()

        info_needed = True
        # TODO - Refactor this so output isn't such an enormous mess. One info field, multiple sub-fields per motif.
        info = '##INFO=<ID=LOCIID,Number=.,Type=String,Description="IDs for loci that variant overlaps.">'
        info += '##INFO=<ID=SAMPSV,Number=.,Type=String,Description="Samples containing the variant.">\n'
        info += ('##INFO=<ID=LOCIVZ,Number=.,Type=String,Description="Z-score for each loci containing the variant.',
                 ' Calculated for each sample containing the variant for each loci.">\n')
        info += ('##INFO=<ID=SAMPTHN,Number=.,Type=Integer,Description="Number of samples in which the variant meets',
                 ' the z-score magnitude threshold.">\n')
        info += ('##INFO=<ID=SAMPSNV,Number=1,Type=Integer,Description="Number of samples containing variant',
                 ' and having loci data.">\n')
        info += ('##INFO=<ID=SAMPSNR,Number=1,Type=Integer,Description="Number of samples containing reference',
                 ' and having loci data.">\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-a", "--activity", dest="activity_file", required=True)
    parser.add_argument("-ov", "--outputvcf", dest="output_vcf", required=True)
    parser.add_argument("-ob", "--outputbed", dest="output_bed", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=0)
    parser.add_argument("-fan", "--filter_a_n", dest="filter_a_n", required=False, default=0)
    # Only print activity if it affects more samples than this number
    # default is -1 so a region with any number of samples affected
    # (having z-scores above threshold) will be output.
    parser.add_argument("-ib", "--include_bed", action="set_true", required=False)
    # Should lines in the activity (bed) output file be include if they don't match a motif?
    # -fa sets this to True
    parser.add_argument("-iv", "--include_vcf", action="set_true", required=False)
    # Should lines in the VCF output files be included if they don't have activity?
    # -fv tag sets this to True

    args = parser.parse_args()

    # Easier to use argument variables
    act_file = args.activity_file
    vcf_out = args.output_vcf
    bed_out = args.output_bed
    inp_file = args.input_file
    th = float(args.threshold)

    filter_bed_num = int(args.filter_a_n)
    include_bed = args.include_bed
    include_vcf = args.include_vcf
    main()
