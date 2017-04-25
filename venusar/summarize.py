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
import time
from math import sqrt

import plotly
import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

from utils import Position, timeString

# YYY-JA 04/24/2017 - Hate making yet another variant of this class. Move to utils.py and make uniform across modules.


class Variant(object):
    """
    Use to process and handle variant records from a VCF more easily. Create from line of VCF file.
    """

    def __init__(self, line):
        self.line_list = line.strip().split("\t")
        self.pos = Position(self.line_list[0], int(self.line_list[1]),
                            (int(self.line_list[1]) + len(self.line_list[3])))
        self.ref_allele = self.line_list[3]
        self.var_allele = self.line_list[4]
        self.iden = self.line_list[2]
        self.orig_line = line.strip()
        self.info_fields = self.line_list[7].split(";")
        (self.common_samples, self.motif_fields, self.exp_fields, self.act_fields,
            self.genes) = self.parse_info_fields()
        self.motif_scores = self.get_motif_scores()

        if self.common_samples is not None:  # Should never evaluate to False.
            self.num_com_samps = len(self.common_samples)
        else:
            self.num_com_samps = 0

    def parse_info_fields(self):
        """
        Get names of samples containing variant and motif INFO fields from a variant record's INFO fields.

        Args:
            self (Variant): Variant object.

        Returns:
            common_samples (dict of tuple): List of samples that had variant called and have loci and expression data.
            motif_fields (list of str): List of INFO fields for variant that contain MOTIF related information.
            exp_fields (list of str): List of INFO fields for variant that contain Expression related information.
            act_fields (list of str): List of INFO fields for variant that contain Activity/Loci related information.
            genes (list of str): List of INFO fields for variant that contain MOTIF related information.
        """
        act_samples = None
        exp_samples = None

        common_samples = []
        motif_fields = []
        exp_fields = []
        act_fields = []
        genes = None

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

            # Break up info fields.
            # YYY-JA 04/25/2017 - could easily be simplified by changing output fields to something more standardized.
            if name.startswith("MOTIF"):
                motif_fields.append(field)
            elif name.startswith("EXP"):
                exp_fields.append(field)
                # Get variant samples with expression data.
                if name == "EXPV":
                    exp_samples = data.split(",")
            elif name.startswith("GENE"):
                genes = data.split(',')
            elif name.startswith("LOCI") or name.startswith("SAMP"):
                act_fields.append(field)
                # Get variant samples with locus activity data.
                if name == "SAMPSTV":
                    act_samples = data.split(",")

        common_samples = compare_samples(exp_samples, act_samples)

        return (common_samples, motif_fields, exp_fields, act_fields, genes)

    def get_motif_scores(self):
        """
        Returns the difference between reference and variant scores for each motif the variant matches
        as a dictionary of {motif_name: (variant - reference log likelihood ratios)}.

        Returns:
            motifs (dict of floats): {motif_name: (variant - reference log likelihood ratios)}
        """

        # List of lists [[field_name1, data1], [field_name2, data2]...]
        motif_fields = [x.split("=") for x in self.motif_fields]
        motifs = {}

        for x in motif_fields:
            # Get motif names and var/ref scores.
            if x[0] == "MOTIFN":
                names = x[1].split(",")
            elif x[0] == "MOTIFV":
                var_scores = x[1].split(",")
            elif x[0] == "MOTIFR":
                ref_scores = x[1].split(",")

        for i in names:
            idx = names.index(i)
            diff = var_scores[idx] - ref_scores[idx]
            motifs[i] = diff

        return motifs

    # XXX-JA 04/25/2017 - Incomplete.

    def get_sample_data(self):
        """
        Parses and returns the gene expression and loci activity z-scores for a given sample.

        Returns:
            motifs (dict of floats): {motif_name: (variant - reference log likelihood ratios)}
        """

        samples = self.common_samples
        act_data = [x.split("=") for x in self.act_fields]

        for x in act_data:
            # Get loci id.
            if x[0] == "LOCIID":
                loci = x[1].split(",")
            # Get each set of z-scores for each loci.
            elif x[0] == "LOCIVZ":
                z_scores = [x.strip("(").strip(")") for x in x[1].split("),(")]  
            

        for i in names:
            idx = names.index(i)
            diff = var_scores[idx] - ref_scores[idx]
            motifs[i] = diff

        return motifs

    # XXX-JA 04/24/2017 - Incomplete. Will create a list of output lines for a given Variant for the
    # summarization report.

    def get_variant_summary(self, d_thresh):
        """
        Creates a list of summary output lines for a given Variant.

        Returns:
            output (list of str): List of output lines for the Variant.
        """


def compare_samples(exp_samples, act_samples):
    """
    Compare gene expression and activity samples with variant and return common samples as a dict
    of format {sample_name: (expression_index, activity_index)} so expression and activity
    z-scores can be found appropriately.

    Args:
        exp_samples (list of str): List of variant samples with expression data.
        act_samples (list of str): List of variant samples with activity data.

    Returns:
        common_samps (dict of tuple): {sample_name: (expression_index, activity_index)}

    """

    common_samps = {}
    samps = list(set(exp_samples) & set(act_samples))  # Get list of common samples.

    # Create new dict for common samples with indexes for sample positions in terms of expression/activity z-scores.
    for x in samps:
        exp_idx = exp_samples.index(x)
        act_idx = act_samples.index(x)
        common_samps[x] = (exp_idx, act_idx)

    return common_samps


def calc_distance(score_array):
    """
    Returns a distance for a motif log-odds ratios, activity z-score, and gene expression z-score.
    """

    return sqrt((score_array[0] ** 2) + (score_array[1] ** 2) + (score_array[2] ** 2))


def plot_distances(score_array, out_prefix):
    """
    Plot distance scores calculated for the motif log-odds ratios, activity z-score, and gene expression
    z-score sets for each variant.

    Args:
        score_array (numpy array) = Array containing all log-odds ratio, activity z-score, and gene expression
            z-score sets for each variant.
        out_prefix (str) = Prefix to use for plot outputs.
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
        d_thresh (float): Distance threshold to be used for more restricted plotting and reporting.
    """

    with open(vcf_file) as f:

        line = f.readline().strip()
        now = time.strftime("%c")

        command = ('##venusaur=<ID=summary,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --expression ' + exp_file + ' --outputvcf ' + out_vcf +
                   ' --threshold ' + str(thresh) + ' --include_vcf ' + str(include_vcf) + '">')

        # Print new info lines at the top of the ##INFO section.
        while line.startswith("##"):
            line = f.readline().strip()

        for line in f:
            current_var = Variant(line)

    print("Complete at: " + timeString() + ".")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_pre", required=True)
    parser.add_argument("-d", "--dthresh", dest="d_thresh", required=False, default=0, type=float)
    args = parser.parse_args()

    main(args.inp_file, args.out_pre, args.d_thresh)
