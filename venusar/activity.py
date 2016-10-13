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
        Default is 0, so all loci a variant overlaps will be reported.
    -fan (int, optional): Set number of samples that must meet z-score threshold for a locus to be reported to
        bed output file. So this number of samples must have the variant and have the locus's activity be significantly
        affected by it. Default is 0, so a locus will be reported if its activity is altered in even one sample above
        the z-score threshold.
    -ib (bool, optional): Should loci that don't contain any variants that significantly affect their activity
        be included in the bed output? False by default, set to True if wanted.
    -iv (bool, optional): Should variants that don't significantly alter a locus's activity be included in the
        vcf output? False by default, set to True if wanted.
"""
import argparse
import time
from statistics import mean, stdev


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


# TODO - Move into a utils.py file and import as appropriate.
class Variant(object):
    """
    Use to process and handle variant records from a VCF more easily. Create from line of VCF file.
    """

    def __init__(self, line, all_sample_names):
        self.line_list = line.strip().split("\t")
        self.pos = Position(self.line_list[0], int(self.line_list[1]), len(self.line_list[3]))
        self.ref_allele = self.line_list[3]
        self.var_allele = self.line_list[4]
        self.iden = self.line_list[2]
        self.orig_line = line.strip()
        self.info_fields = self.line_list[7].split(";")
        self.var_samples, self.motif_fields = self.parse_info_fields()
        self.ref_samples = [x for x in all_sample_names if x not in self.var_samples]
        self.loci = []

        if self.var_samples is not None:  # Should never evaluate to False.
            self.num_var_samps = len(self.samples)
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
            field_info = field.split("=")
            name, data = (field_info[0], field_info[1])

            # TODO - Write method that parses header to determine # samples with variant rather than this lazy method.
            if name is "set":
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
        info.insert(0, "SAMPSR=" + ",".join([x.ref_samples[self] for x in self.loci[0]]))
        info.insert(0, "SAMPSV=" + ",".join([x.var_samples[self] for x in self.loci[0]]))
        info.insert(0, "SAMPSNR=" + ",".join([x.num_valid_ref[self] for x in self.loci[0]]))
        info.insert(0, "SAMPSNV=" + ",".join([x.num_valid_var[self] for x in self.loci[0]]))

        # Use lists to maintain order in output so that LOCIID, LOCIVZ, SAMPTHN fields can all be matched up.
        z_scores = []
        pass_thresh = []
        loci_idens = []
        for item in self.loci:
            loci_idens.append(item.iden)
            pass_thresh.append(item.num_pass_thresh[self])
            tmp = "(" + ",".join([x for x in item.z_scores[self]]) + ")"
            z_scores.append(tmp)
        info.insert(0, "SAMPTHN=" + ",".join(pass_thresh))
        info.insert(0, "LOCIVZ=" + ",".join(z_scores))
        info.insert(0, "LOCIID=" + ",".join(loci_idens))

        if include_vcf and all([x == 0 for x in pass_thresh]):
            return None
        else:
            self.info_fields = info
            self.line_list[7] = ";".join(self.info_fields)
            output = "\t".join(self.line_list)
            return output


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
        self.pos = Position(line_list[0], int(line_list[1]), len(line_list[2]))
        self.orig_line = line.strip()
        self.iden = str(line[3])
        self.data = [float(x) for x in line[4:]]
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
        Add Variant object variant to list of variants that overlap the Locus.
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
        self.num_valid_ref = len(ref_ind)
        self.num_valid_var = len(var_ind)

        if len(ref_ind) > 1:  # If all samples (or all but 1) have the variant, can't calc z-score, return 'NA'.
            for entry in ref_ind:
                scores = self.ref_scores[variant]
                scores.append(self.data[int(entry)])
                self.ref_scores[variant] = scores
            for entry in var_ind:
                scores = self.var_scores[variant]
                scores.append(self.data[int(entry)])
                self.var_scores[variant] = scores

            ref_mean = mean(self.ref_scores)
            ref_std = stdev(self.ref_scores)

            if ref_std != 0:  # If only one sample has ref, will have no variance. Should never happen.
                for i in var_ind:
                    vals = self.z_scores[variant]
                    vals.append("NA")
                    self.z_scores[variant] = vals
            else:
                score = [((x - ref_mean) / ref_std) for x in self.var_scores]
                if abs(score) >= thresh:  # Check number of variant samples that passed the threshold.
                    passed = self.num_pass_thresh[variant]
                    passed += 1
                    self.num_pass_thresh[variant]
                vals = self.z_scores[variant]
                vals.append(score)
                self.z_scores[variant] = vals
        else:
            for i in var_ind:
                vals = self.z_scores[variant]
                vals.append("NA")
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

            out_line = [chrom, start, end, iden]
            motifs_out = ";".join(item.motif_fields)
            out_line.append(item.pos.chrom + ":" + item.pos.start + "_" + item.ref_allele + ">" + item.var_allele)
            val_var_samps = self.var_samples[item]
            num_val_var = len(self.var_samples[item])
            val_ref_samps = ",".join(self.ref_samples[item])
            num_val_ref = len(self.ref_samples[item])
            all_var_samps = ",".join(item.var_samples)
            num_all_var_samps = len(item.var_samples)
            all_ref_samps = ",".join(item.ref_samples)
            num_all_ref_samps = len(item.var_samples)

            # Handle z_scores.
            scores_out = []
            scores = self.z_scores[item]
            for x in range(len(val_var_samps)):
                samp = val_var_samps[x]
                score = scores[x]
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


def main(vcf_file, act_file, out_vcf, out_bed, thresh=0, filter_num=0, include_bed=False, include_vcf=False):
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
    """
    print("Parsing activity data file.")
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
        info += ('##INFO=<ID=LOCIVZ,Number=.,Type=String,Description="Z-score for each loci containing the variant.'
                 ' Calculated for each sample containing the variant for each loci.">\n')
        info += ('##INFO=<ID=SAMPTHN,Number=.,Type=Integer,Description="Number of samples in which the variant meets'
                 'the z-score magnitude threshold.">\n')
        info += ('##INFO=<ID=SAMPSNV,Number=1,Type=Integer,Description="Number of samples containing variant'
                 ' and having loci data.">\n')
        info += ('##INFO=<ID=SAMPSNR,Number=1,Type=Integer,Description="Number of samples containing reference'
                 ' and having loci data.">')
        command = ('##venusaur=<ID=activity,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --activity ' + act_file + ' --outputvcf ' + out_vcf + ' --outputbed ' + out_bed +
                   ' --threshold ' + str(thresh) + ' --filter_act_num ' + str(filter_num) + ' --include_bed ' +
                   include_bed + ' --include_vcf ' + include_vcf + '">')

        # Print new info lines at the top of the ##INFO section.
        while line.startswith("##"):
            if info_needed and line.startswith("##INFO"):
                print(command, file=output_vcf)
                print(command, file=output_bed)
                print(info, file=output_vcf)
                info_needed = False
            print(line, file=output_vcf)
            line = f.readline().strip()

        vcf_samples = get_vcf_samples(line)  # Parse VCF sample header line to get samples present in file.

        print("Comparing samples in VCF file and activity file to find commonalities.\n")
        print("VCF samples: ", *vcf_samples, sep="\n")
        print("Activity samples: ", *list(act_samps.keys()), sep="\n")

        common_samps, valid_act_samps = compare_samples(act_samps, vcf_samples)  # Get common samples b/twn the two.
        print("Common samples: ", *common_samps, sep="\n")
        print("Processing variants. This may take some time.")
        # TODO - Progress bar might actually be a decent addition.

        for line in f:
            current_var = Variant(line, vcf_samples)
            loci_ovlp_var = []

            for item in act_data:
                if current_var.pos.chrom == item.pos.chrom and current_var.pos.start > item.pos.end:
                    break
                elif current_var.pos.chrom != item.pos.chrom:
                    continue
                elif current_var.pos.overlaps(item.pos):
                    loci_ovlp_var.append(item)

            # If variant overlaps no loci, print to output only if include_vcf option used.
            if not loci_ovlp_var:
                if include_vcf:
                    print(line, file=output_vcf)
                    continue
                else:
                    continue

            # Get activity data indices for both samples with variant and without.
            var_act_indices = [valid_act_samps[x] for x in current_var.samples if x in valid_act_samps]
            ref_act_indices = [valid_act_samps[x] for x in valid_act_samps if x not in current_var.samples]

            # Calculate z-scores.
            for x, loc in enumerate(loci_ovlp_var):
                var_samples = [x for x in current_var.samples if x in valid_act_samps]
                ref_samples = [x for x in valid_act_samps if x not in current_var.samples]
                loc.add_variant(current_var, var_samples, ref_samples)  # Add Variant to Locus object.
                loc.calc_z_score(ref_act_indices, var_act_indices, current_var, thresh)
                current_var.loci.append(loc)  # Add Locus object to given Variant.
                loci_ovlp_var[x] = loc
                loci_out.append(loc)  # These will be used for eventual BED output.

            vcf_out_line = current_var.get_variant_output(include_vcf)

            if vcf_out_line is not None:
                print(vcf_out_line, file=output_vcf)

    print("Filtering loci and creating BED output.")
    print("CHR", "START", "END", "ID", "VARIANT", "Z_SCORES", "NUM_PASS_THRESH", "COMMON_VAR_SAMPS",
          "NUM_COMMON_VAR_SAMPS", "COMMON_REF_SAMPS", "NUM_COMMON_REF_SAMPS", "ALL_VAR_SAMPS", "ALL_REF_SAMPS",
          "MOTIF_INFO", sep="\t", file=output_bed)
    for item in loci_out:
        bed_out_line = item.get_locus_output(include_bed, filter_num)

        if bed_out_line is not None:
            print(*bed_out_line, sep="\n", file=output_bed)

    print("Complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-a", "--activity", dest="activity_file", required=True)
    parser.add_argument("-ov", "--outputvcf", dest="output_vcf", required=True)
    parser.add_argument("-ob", "--outputbed", dest="output_bed", required=True)
    parser.add_argument("-th", "--threshold", dest="threshold", required=False, default=0)
    parser.add_argument("-fan", "--filter_act_num", dest="filter_a_n", required=False, default=0)
    parser.add_argument("-ib", "--include_bed", action="set_true", required=False)
    parser.add_argument("-iv", "--include_vcf", action="set_true", required=False)

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
