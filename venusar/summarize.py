#!/usr/bin/env python
"""
For a VEGAS'd VCF file, calculate a distance score for each motif, ChIP z-score, and gene expression z-score
set for each variant. These values will be plotted in 3 dimensional space. Additionally, simple bed-like
files are provided as output. One contains the scores for all motif, ChIP z-score, and gene expression z-score
sets for all variants. An optional second contains only those sets that meet the distance score threshold as
defined by the user. A third will report the sets for the top 100 distance scores.

Usage: summarize.py -i <input.vcf> -o <output> [OPTIONS]

Args:
    -i (str): Path to sorted variant file to process.
    -o (str): Prefix for output files.
    -d (float, optional): Distance magnitude threshold that must be met for variants/genes to be reported to output.
        Default is 0, so all variant-sample activity-gene set distances will be reported.
"""
import argparse

import plotly
import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

from sequence import read_line2sample_list
from utils import Position

# YYY-JA 04/24/2017 - Hate making yet another variant of this class. Move to utils.py and make uniform across modules.


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
        (self.var_samples, self.motif_fields, self.exp_fields, self.act_fields,
            self.gene_fields) = self.parse_info_fields()
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
        exp_fields = []
        gene_fields = []

        for field in self.info_fields:
            if field != "INDEL":  # Take care of INDEL flag.
                field_info = field.split("=")

                # YYY-JA - This is a hack around a bug that's messing up the MOTIFN field in tf_expression.py.
                # Go back and actually figure out why the MOTIFN field is getting split up sometimes.
                try:
                    name, data = (field_info[0], field_info[1])
                except:
                    name, data = "BROKEN", None
            else:
                name, data = "INDEL", None

            # YYY-JA - Write method to parse header to determine # samples with variant rather than this lazy method.
            if name == "SAMPSTV":
                samples = data.split(",")
            elif name.startswith("MOTIF"):
                motif_fields.append(field)
            elif name.startswith("EXP"):
                exp_fields.append(field)
            elif name.startswith("GENE"):
                gene_fields.append(field)

        return (samples, motif_fields, exp_fields, gene_fields)


    # XXX-JA 04/24/2017 - Incomplete. Will create a list of output lines for a given Variant for the
    # summarization report.
    def get_variant_summary(self, d_thresh):
        """
        Creates a list of summary output lines for a given Variant
        """


def plot_distances(score_array, out_prefix, distance_threshold=0):
    """
    Plot distance scores calculated for the motif log-odds ratios, activity z-score, and gene ex-ression
    z-score sets for each variant.

    Args:
        score_array (numpy array) = Array containing all log-odds ratio, activity z-score, and gene expression
            z-score sets for each variant.
        out_prefix (str) = Prefix to use for plot outputs.
        distance_threshold (float) = Value used to filter distances from being plotted. Used as a magnitude value,
            anything less will be excluded from the plot and second optional output file.
    """

    x, y, z = np.random.multivariate_normal(np.array([0, 0, 0]), np.eye(3), 400).transpose()

    trace1 = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=4,                 # Using gene expression as color scale values for now.
            color=z,                # set color to an array/list of desired values
            colorscale='Viridis',   # choose a colorscale
            opacity=0.8
        )
    )

    data = [trace1]
    layout = go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        )
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=out_prefix + '.html')
    py.image.save_as(fig, filename=out_prefix + '.png')


def main(vcf_file, out_prefix, d_thresh):
    """

    Args:
        vcf_file (str): Path to sorted variant file to process.
        out_prefix (str): Prefix to be used for output files.
        d_thresh (float):  
    """

    with open(vcf_file) as f:

        line = f.readline().strip()
        now = time.strftime("%c")

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

        vcf_samples = read_line2sample_list(line)  # Parse VCF sample header line to get samples present in file.

        for line in f:
            current_var = Variant(line, vcf_samples)
            loci_ovlp_var = []

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

    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_pre", required=True)
    parser.add_argument("-d", "--dthresh", dest="d_thresh", required=False, default=0, type=float)
    args = parser.parse_args()

    main(args.inp_file, args.out_pre, args.d_thresh)
