#!/usr/bin/env python
"""
For a given motif annotated VCF file (already run through motifs.py) and a bed-like file for loci of interest and some
value for each loci for each sample, find loci that overlap a variant and compare the value of samples with the variant
to those without the variant. Report robut z-scores for each loci overlapped in an output VCF and report the variants
for each loci in a bed-like, loci-centric output file as well.

Usage: activity.py -i <input.vcf> -a <activity.bed> -ov <output.vcf> -ob <output.bed> [OPTIONS]

Args:
    -i (str): Path to sorted variant file to process.
    -a (str): Path to activity 'bed' file.
    -ov (str): Path to VCF output file to be created.
    -ob (str): Path to loci output file to be created.
    -th (float, optional): Z-score magnitude threshold that must be met for variants/loci to be reported to output.
        Default is 0, so all loci a variant overlaps will be reported.
    -fan (int, optional): Set number of samples that must meet z-score threshold for a locus to be reported to
        bed output file. So this number of samples must have the variant and have the locus's activity be significantly
        affected by it. Default is 0, so a locus will be reported if its activity is altered in even one sample above
        the robust z-score threshold.
    -ib (bool, optional): Should loci that don't contain any variants that significantly affect their activity
        be included in the bed output? False by default, set to True if wanted.
    -iv (bool, optional): Should variants that don't significantly alter a locus's activity be included in the
        vcf output? False by default, set to True if wanted.
"""

from __future__ import print_function    # so Ninja IDE will stop complaining & show symbols
import argparse
import time
from statistics import median
from sequence import read_line2sample_list
from utils import Position, timeString


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
            if name == "set":
                samples = data.split("-")
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
        info.insert(0, "SAMPSTV=" + ",".join(self.var_samples))
        info.insert(0, "SAMPSR=" + ",".join([x.ref_samples[self] for x in self.loci][0]))
        info.insert(0, "SAMPSV=" + ",".join([x.var_samples[self] for x in self.loci][0]))

        # TODO - Check and make sure next two lines are functioning properly.
        info.insert(0, "SAMPSNR=" + ",".join([str(x.num_valid_ref[self]) for x in self.loci]))
        info.insert(0, "SAMPSNV=" + ",".join([str(x.num_valid_var[self]) for x in self.loci]))

        # Use lists to maintain order in output so that LOCIID, LOCIVZ, SAMPTHN fields can all be matched up.
        z_scores = []
        pass_thresh = []
        loci_idens = []
        for item in self.loci:
            loci_idens.append(item.iden)
            pass_thresh.append(item.num_pass_thresh[self])
            tmp = "(" + ",".join([str(round(x, 4)) for x in item.z_scores[self][0]]) + ")"
            z_scores.append(tmp)
        info.insert(0, "SAMPTHN=" + ",".join([str(x) for x in pass_thresh]))
        info.insert(0, "LOCIVZ=" + ",".join(z_scores))
        info.insert(0, "LOCIID=" + ",".join(loci_idens))

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
        """
        Calculate a robust z-score for the given locus and variant.

        This uses the median absolute deviation (MAD):
        https://en.wikipedia.org/wiki/Median_absolute_deviation
        """
        self.num_valid_ref[variant] = len(ref_ind)
        self.num_valid_var[variant] = len(var_ind)

        for entry in ref_ind:
            scores = self.ref_scores[variant]
            scores.append(self.data[int(entry)])
            self.ref_scores[variant] = scores
        for entry in var_ind:
            scores = self.var_scores[variant]
            scores.append(self.data[int(entry)])
            self.var_scores[variant] = scores

        # MAD calculation.
        all_scores = self.ref_scores[variant] + self.var_scores[variant]
        med = median(all_scores)
        abs_score = [abs(x - med) for x in all_scores]
        mad = median(abs_score) * 1.4826
        # 1.4826 is a constant that assumes a normal distribution to use the MAD as a consistent estimator
        # of standard deviation.

        robust_z_scores = [((x - med) / mad) for x in self.var_scores[variant]]

        for item in robust_z_scores:
            if abs(item) >= thresh:  # Check number of variant samples that passed the threshold.
                passed = self.num_pass_thresh[variant]
                passed += 1
                self.num_pass_thresh[variant] = passed
            vals = self.z_scores[variant]
            vals.append(robust_z_scores)
            self.z_scores[variant] = vals

        return

    def get_locus_output(self, include_bed=False, filter_num=0):
        """
        Create list of output lines for given Locus object.

        Args:
            include_bed (bool): True if variants that don't pass the z-score threshold for any Locus should excluded
                from output. False if they should be included.
            filter_num (int): Number of samples the must meet z-score threshold for Variant to be included in output.

        Returns:
            output (list of str): List of lines for Locus in appropriate BED-like format.
                or
            None: If include_bed is True and Locus doesn't contain any Variants that hit z-score threshold for any
                sample.
        """
        output = []
        chrom = self.pos.chrom
        start = self.pos.start
        end = self.pos.end
        iden = self.iden

        # Get info for each Variant that is necessary for output.
        for item in self.variants:
            num_meet_thresh = int(self.num_pass_thresh[item])
            if num_meet_thresh < filter_num and include_bed is False:
                continue  # Just go to next Variant if required number of samples didn't meet z-score threshold.

            out_line = [chrom, str(start), str(end), str(iden)]
            motifs_out = ";".join(item.motif_fields)
            out_line.append(item.pos.chrom + ":" + str(item.pos.start) + "_" + item.ref_allele + ">" + item.var_allele)
            val_var_samps = self.var_samples[item]
            num_val_var = len(self.var_samples[item])
            val_ref_samps = ",".join(self.ref_samples[item])
            num_val_ref = len(self.ref_samples[item])
            all_var_samps = ",".join(item.var_samples)
            num_all_var_samps = len(item.var_samples)
            all_ref_samps = ",".join(item.ref_samples)
            num_all_ref_samps = len(item.ref_samples)

            # Handle z_scores.
            scores_out = []
            scores = self.z_scores[item]
            for x in range(len(val_var_samps)):
                samp = val_var_samps[x]
                score = round(scores[x][x], 4)
                scores_out.append(samp + "=" + str(score))

            out_line.append(",".join(scores_out))
            out_line.append(str(num_meet_thresh))
            out_line.append(",".join(val_var_samps))
            out_line.append(str(num_val_var))
            out_line.append(val_ref_samps)
            out_line.append(str(num_val_ref))
            out_line.append(all_var_samps)
            out_line.append(str(num_all_var_samps))
            out_line.append(all_ref_samps)
            out_line.append(str(num_all_ref_samps))
            out_line.append(motifs_out)

            output.append("\t".join(out_line))

        if output:  # If no variants are in the output list, just return None.
            return output
        else:
            return None


def get_activity_samples(header_line):
    """
    Parse header of activity file to return sample names and column indices.

    Args:
        header_line (str): Header line from activity file.

    Returns:
        act_samples (dict): Dictionary of {sample_name (str): sample_data_index (int)}.
            sample_data_index is index for data in sample list, not the line as a whole.
            e.g.: [samp1, samp2, samp3] & [20, 10, 5] for data values, then {'samp1': 0}.
    """
    line_list = header_line.strip().split("\t")
    samples = line_list[4:]
    act_samples = {}

    for item in samples:
        samp_idx = samples.index(item)
        sample = item.split(".")[0]
        act_samples[sample] = samp_idx

    return act_samples


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
        header = f.readline().strip()
        act_samples = get_activity_samples(header)  # Get sample names/indices.
        act_data = []

        for line in f:
            record = Locus(line)
            act_data.append(record)

        return (act_samples, act_data)


def reduce_activity_names(act_samps, split_string="_"):
    """
    Return only unique part of names in the passed list set.

    Code assumes either start or end is unique based on unique set size

    Args:
        act_samps: list of strings for the activity samples
            or dictionary with keys being strings for activity samples
        split_string: String to split individual act_samps on
            assumes single split is relevant
            default = "_"
    Returns:
        act_samps modified to not include the non-unique part of the input strings
        returns same type as the input act_samps, list or dict
    """

    split_one = []
    split_two = []

    for sample in act_samps:
        # 1. split on split_string
        splitList = sample.split(split_string)
        # 2. put first split and all remaining splits in 2 arrays
        split_one.append(splitList[0])
        if len(splitList) == 1:
            # because otherwise it adds an empty list item and breaks below
            split_two.append("")
        else:
            split_two.append(split_string.join(splitList[1:]))

    # 3. determine the unique set size (just making it a set makes them unique
    s1 = set(split_one)
    s2 = set(split_two)

    if len(s1) > len(s2):
        # s2 is the non-unique part; ie s1 is unique
        act_samps_temp = list(s1)
    else:
        # s1 is the non-unique part; ie s2 is unique
        act_samps_temp = list(s2)

    if type(act_samps) is list:
        # do nothing just return
        return (act_samps_temp)
    elif type(act_samps) is dict:
        # must rebuild the dictionary
        act_samps_rebuild = {}
        ind = -1
        for sample in act_samps:
            ind = ind + 1
            act_samps_rebuild[act_samps_temp[ind]] = act_samps[sample]
        return (act_samps_rebuild)


def main(vcf_file, act_file, out_vcf, out_bed, thresh=0, filter_num=0, include_bed=False, include_vcf=False,
         drop_act_=1):
    """
    Compare activity of loci for samples harboring a variant within a given locus to those samples that do not.

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
        drop_act_ (integer, optional): If > 0 then break activity items on _,
            return only unique part of name.
            code assumes either start or end is unique based on unique set size
            once dropped reruns comparison to the vcf samples
            if 1: only runs if prior vcf comparison results in no overlap
            if 2: runs no matter what
    """
    print("Parsing activity data file: " + timeString() + ".")
    act_samps, act_data = parse_activity_file(act_file)

    output_vcf = open(out_vcf, "w")
    output_bed = open(out_bed, "w")

    loci_out = []  # Use to hold all Locus objects that overlap a Variant.

    with open(vcf_file) as f:

        # Add new INFO lines.
        line = f.readline().strip()
        now = time.strftime("%c")
        info_needed = True
        # TODO - Refactor this so output isn't such an enormous mess. One info field, multiple sub-fields per motif.
        # TODO - Add sample names for those that pass threshold.
        info = '##INFO=<ID=LOCIID,Number=.,Type=String,Description="IDs for loci that variant overlaps.">\n'
        info += '##INFO=<ID=SAMPSTV,Number=.,Type=String,Description="All samples with the variant allele.">\n'
        info += ('##INFO=<ID=SAMPSR,Number=.,Type=String,Description="Samples with the reference allele and loci data'
                 '.">\n')
        info += ('##INFO=<ID=SAMPSV,Number=.,Type=String,Description="Samples with the variant allele and loci data.'
                 '">\n')
        info += ('##INFO=<ID=LOCIVZ,Number=.,Type=String,Description="Robust z-score for each loci '
                 'containing the variant. Calculated for each sample containing the variant for each loci.">\n')
        info += ('##INFO=<ID=SAMPTHN,Number=.,Type=Integer,Description="Number of samples in which the variant meets'
                 'the z-score magnitude threshold.">\n')
        info += ('##INFO=<ID=SAMPSNV,Number=1,Type=Integer,Description="Number of samples containing variant'
                 ' and having loci data.">\n')
        info += ('##INFO=<ID=SAMPSNR,Number=1,Type=Integer,Description="Number of samples containing reference'
                 ' and having loci data.">')
        command = ('##venusaur=<ID=activity,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --activity ' + act_file + ' --outputvcf ' + out_vcf + ' --outputbed ' + out_bed +
                   ' --threshold ' + str(thresh) + ' --filter_act_num ' + str(filter_num) + ' --include_bed ' +
                   str(include_bed) + ' --include_vcf ' + str(include_vcf) + '">')

        # Print new info lines at the top of the ##INFO section.
        while line.startswith("##"):
            if info_needed and line.startswith("##INFO"):
                print(command, file=output_vcf)
                print(command, file=output_bed)
                print(info, file=output_vcf)
                info_needed = False
            print(line, file=output_vcf)
            line = f.readline().strip()

        vcf_samples = read_line2sample_list(line)  # Parse VCF sample header line to get samples present in file.

        print(line, file=output_vcf)

        print("Comparing samples in VCF file and activity file to find commonalities.\n")
        print("VCF samples: ", *vcf_samples, end="\n\n")
        print("Activity samples: ", *list(act_samps.keys()), end="\n\n")

        common_samps, valid_act_samps = compare_samples(act_samps, vcf_samples)  # Get common samples b/twn the two.
        print("Common samples: ", *common_samps, end="\n\n")
        if drop_act_ > 0:
            if drop_act_ == 1 and len(common_samps) == 0:
                redo_compare = True
                act_samps = reduce_activity_names(act_samps)
            elif drop_act_ == 2:
                redo_compare = True
                # merge old and new samps to match when compare_samples is run below
                # if they were just lists the following would work but they are not
                # act_samps = list(set(reduce_activity_names(act_samps)) | set(list(act_samps)))
                extend_dict = reduce_activity_names(act_samps)
                for extdictkey in extend_dict:
                    act_samps[extdictkey] = extend_dict[extdictkey]
            else:
                redo_compare = False
            if redo_compare:
                # Get common samples b/twn the two input sets: vcf and activity.
                common_samps, valid_act_samps = compare_samples(act_samps, vcf_samples)
                print("Updated Common samples: ", *common_samps, end="\n\n")

        print("Processing variants. This may take some time.")
        # TODO - Progress bar might actually be a decent addition.

        for line in f:
            current_var = Variant(line, vcf_samples)
            loci_ovlp_var = []

            # Check if any of the variant samples actually have activity data as well, skip if not.
            for x in current_var.var_samples:
                if x in common_samps:
                    for item in act_data:
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

            # Get activity data indices for both samples with variant and without.
            var_act_indices = [valid_act_samps[x] for x in current_var.var_samples if x in valid_act_samps]
            ref_act_indices = [valid_act_samps[x] for x in valid_act_samps if x not in current_var.var_samples]

            # Calculate z-scores.
            for x, loc in enumerate(loci_ovlp_var):
                var_samples = [x for x in current_var.var_samples if x in valid_act_samps]
                ref_samples = [x for x in valid_act_samps if x not in current_var.var_samples]
                loc.add_variant(current_var, var_samples, ref_samples)  # Add Variant to Locus object.
                loc.calc_z_score(ref_act_indices, var_act_indices, current_var, thresh)
                current_var.loci.append(loc)  # Add Locus object to given Variant.
                loci_ovlp_var[x] = loc
                if loc not in loci_out:
                    loci_out.append(loc)  # These will be used for eventual BED output.

            vcf_out_line = current_var.get_variant_output(include_vcf)

            if vcf_out_line is not None:
                print(vcf_out_line, file=output_vcf)
            elif include_vcf:
                print(line.strip(), file=output_vcf)

    print("Filtering loci and creating BED output.")
    print("CHR", "START", "END", "ID", "VARIANT", "Z_SCORES", "NUM_PASS_THRESH", "COMMON_VAR_SAMPS",
          "NUM_COMMON_VAR_SAMPS", "COMMON_REF_SAMPS", "NUM_COMMON_REF_SAMPS", "ALL_VAR_SAMPS", "NUM_ALL_VAR_SAMPS",
          "ALL_REF_SAMPS", "NUM_COMMON_REF_SAMPS"
          "MOTIF_INFO", sep="\t", file=output_bed)
    for item in loci_out:

        bed_out_line = item.get_locus_output(include_bed, filter_num)

        if bed_out_line is not None:
            print(*bed_out_line, sep="\n", file=output_bed)

    print("Complete: " + timeString() + ".")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-a", "--activity", dest="activity_file", required=True)
    parser.add_argument("-ov", "--outputvcf", dest="output_vcf", required=True)
    parser.add_argument("-ob", "--outputbed", dest="output_bed", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=0)
    parser.add_argument("-fan", "--filter_act_num", dest="filter_a_n", required=False, default=0)
    parser.add_argument("-ib", "--include_bed", action="store_true", required=False)
    parser.add_argument("-iv", "--include_vcf", action="store_true", required=False)

    args = parser.parse_args()

    inp_file = args.input_file
    act_file = args.activity_file
    vcf_out = args.output_vcf
    bed_out = args.output_bed
    th = float(args.threshold)
    filter_bed_num = int(args.filter_a_n)
    include_bed = args.include_bed
    include_vcf = args.include_vcf

    main(inp_file, act_file, vcf_out, bed_out, th, filter_bed_num, include_bed, include_vcf)
