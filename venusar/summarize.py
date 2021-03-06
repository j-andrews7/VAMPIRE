#!/usr/bin/env python
"""
For a VEGAS'd VCF file, calculate a distance score for each motif delta score, ChIP z-score, and gene expression
z-score set for each variant. These values will be plotted in 3 dimensional space. Additionally, simple bed-like files
are provided as output. One contains the scores for all motif, ChIP z-score, and gene expression z-score sets for all
variants. An optional second contains only those sets that meet the distance score threshold as defined by the user.
A third will report the sets for the top 100 distance scores.

Usage: summarize.py -i <input.vcf> -o <output> [OPTIONS]

Args:
    -i (str): Path to sorted variant file to process.
    -o (str): Prefix for output files.
    -d (float, optional): Distance magnitude threshold that must be met for variants/genes to be reported to output.
        Default is 0, so all variant-sample activity-gene set distances will be reported.
"""
import argparse
import time
import pandas as pd
from math import sqrt

import plotly
import plotly.graph_objs as go

import numpy as np

from utils import Position, timeString

# YYY-JA 04/24/2017 - Hate making yet another variant (Ha) of this class.
# Move to utils.py and make uniform across modules.


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
        self.motif_scores = self.get_motif_scores()  # Calculate delta motif scores and place into a dict.
        self.sample_data = self.get_sample_data()  # Parse all combos of gene expression and loci activity data.

        self.output = self.get_variant_summary()  # Get output lines as a list of lists.

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
                if name == "SAMPSV":
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
            diff = float(var_scores[idx]) - float(ref_scores[idx])
            motifs[i] = diff

        return motifs

    def get_sample_data(self):
        """
        Parses and returns the gene expression and loci activity z-scores for variant samples.

        Returns:
            sample_data (list of lists of str): [[sample, loci1, sample act_z_score, gene1, sample exp_z_score],
                                                 [sample, loci1, sample act_z_score, gene2, sample exp_z_score]...]
        """

        samples = self.common_samples
        act_data = [x.split("=") for x in self.act_fields]
        exp_data = [x.split("=") for x in self.exp_fields]
        genes = self.genes

        sample_data = []

        for x in act_data:
            # Get loci id.
            if x[0] == "LOCIID":
                loci = x[1].split(",")
            # Get each set of z-scores for each loci.
            elif x[0] == "LOCIVZ":
                act_z_scores = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]

        # Create dict for loci - {loci_id: [act_z_scores]}
        loci_data = {k: v for k, v in zip(loci, act_z_scores)}

        for x in exp_data:
            # Get each set of z-scores for each gene.
            if x[0] == "EXPVZ":
                gene_z_scores = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]

        # Create dict for genes - {gene: [gene_z_scores]}
        gene_data = {k: v for k, v in zip(genes, gene_z_scores)}

        # Create list of lists containing all sets of relevant data combinations for all samples with a variant.
        # [[sample, loci1, sample act_z_score, gene1, sample exp_z_score],
        #  [sample, loci1, sample act_z_score, gene2, sample exp_z_score]...]
        for i in samples:  # First iterate through all samples and get their expression and activity indices.
            e_idx = samples[i][0]
            a_idx = samples[i][1]

            for l in loci_data:  # Iterate through each locus for given sample.
                # print(loci_data[l], str(a_idx), sep="\t")
                samp_act_z = loci_data[l][a_idx]

                for g in gene_data:  # And now each gene for each locus.
                    samp_exp_z = gene_data[g][e_idx]

                    sample_data.append([i, l, samp_act_z, g, samp_exp_z])

        return sample_data

    def get_variant_summary(self):
        """
        Creates a list of summary output fields for a given Variant as well as its distance metrics for plotting.

        Returns:
            output_fields (list of lists of str): List of lists of output fields for the Variant.
        """
        var_info = [self.pos.chrom, self.pos.start, self.ref_allele, self.var_allele]
        motif_info = self.motif_scores
        sample_info = self.sample_data

        output_fields = []

        for m in motif_info:
            m_score = motif_info[m]  # Motif delta score for var vs ref.

            for s in sample_info:
                dist_metrics = [(float(m_score)), float(s[2]), float(s[4])]  # To be used for plotting later.
                dist_score = calc_distance(dist_metrics)
                # Create complete list of output fields.
                output_fields.append(var_info + [m] + [m_score] + [s[0]] + [s[1]] + [s[2]] + [s[3]] +
                                     [s[4]] + [dist_score])

        return (output_fields)


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


def plot_distances(df, out_prefix):
    """
    Plot distance scores calculated for the delta var/ref motif log-odds ratios, activity z-score, and gene expression
    z-score sets for each variant.

    Args:
        df (pandas Dataframe) = Dataframe containing all variant, distance metric, and sample info. One record per row.
        out_prefix (str) = Prefix to use for plot outputs.
    """

    # Take top 30k hits only, plotly handle really handle more for 3D plots.
    if len(df) > 30000:
        df = df.head(30000)

    info = list(zip(df.SAMPLE, df.MOTIF, df.GENE))
    info_list = ["Sample: " + x[0] + ", Motif: " + x[1] + ", Gene: " + x[2] for x in info]

    trace1 = go.Scatter3d(
        name="Distances",
        x=df['VAR-REF_SCORE'],
        y=df.ACT_ZSCORE,
        z=df.EXP_ZSCORE,
        hoverinfo="x+y+z+text",
        text=info_list,
        mode='markers',
        marker=dict(
            size=4,                 # Using gene expression as color scale values for now.
            color=df.SCALED_DISTANCE,                # set color to an array/list of desired values
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
            t=100
        ),
        title='Distances from Average for Individual Variant Events',
        scene=dict(
            xaxis=dict(
                title='Var/Ref Motif Log-Odds Ratio Difference',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            ),
            yaxis=dict(
                title='Enhancer Activity z-Score',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            ),
            zaxis=dict(
                title='Gene Expression z-Score',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    # py.image.save_as(fig, filename=out_prefix + '.png')
    plotly.offline.plot(fig, filename=out_prefix + '.html', auto_open=False, image="png", image_filename=out_prefix)


def scale_and_frame(all_output):
    """
    Return dataframe after adding scaled distance score to each row. Scaling done by dividing with max value of each
    metric, summation, and taking square root. Make things much easier to plot and handle.

    Args:
        all_output (list of lists): Each list contains a set of motif, expression, activity data.
            [sample, loci1, sample act_z_score, gene1, sample exp_z_score, distance]

    Returns:
        df (pandas Dataframe): Each row contains a set of motif, expression, activity data.
            [sample, loci1, sample act_z_score, gene1, sample exp_z_score, distance, scaled_distance]
    """

    df = pd.DataFrame(all_output, columns=['CHR', 'POS', 'REF', 'ALT', 'MOTIF', 'VAR-REF_SCORE', 'SAMPLE', 'LOCIID',
                                           'ACT_ZSCORE', 'GENE', 'EXP_ZSCORE', 'DISTANCE'])
    df[['VAR-REF_SCORE', 'ACT_ZSCORE', 'EXP_ZSCORE', 'DISTANCE']] = df[['VAR-REF_SCORE', 'ACT_ZSCORE', 'EXP_ZSCORE',
                                                                        'DISTANCE']].apply(pd.to_numeric)
    print(df.head(10))

    # Get scaling factors.
    motif_score_scale = max(df['VAR-REF_SCORE'] ** 2)
    act_score_scale = max(df.ACT_ZSCORE ** 2)
    gene_score_scale = max(df.EXP_ZSCORE ** 2)

    # Scale the distance and create new column.
    df['SCALED_DISTANCE'] = np.sqrt(((df['VAR-REF_SCORE'] ** 2) / motif_score_scale) +
                                    ((df.ACT_ZSCORE ** 2) / act_score_scale) +
                                    ((df.EXP_ZSCORE ** 2) / gene_score_scale))

    # Sort by distance.
    df.sort_values('SCALED_DISTANCE', ascending=False, inplace=True)

    return df


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

        command = ('##venusar=<ID=summary,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --output ' + out_prefix + ' --dthresh ' + str(d_thresh) + '">')
        print(command)

        full_out_file = open(out_prefix + "_full.txt", "w")  # Full output.
        top_out_file = open(out_prefix + "_top100.txt", "w")  # Top 100 hits.

        if d_thresh != 0:
            rest_out_file = open(out_prefix + "_restricted.txt", "w")  # Restricted by distance threshold.

        all_output = []

        # Print new info lines at the top of the ##INFO section.
        while line.startswith("#"):
            line = f.readline().strip()

        print("Reading and processing input file.")
        for line in f:
            current_var = Variant(line)  # Most processing happening here.

            # Piece out full hits, restricted hits, top hits for everything.
            full_var_output = current_var.output

            for x in full_var_output:
                all_output.append(x)

        # Scale distance metrics.
        print("Calculating distance metrics.")
        scaled_df = scale_and_frame(all_output)

        print("Creating output files and plots.")
        scaled_df.to_csv(rest_out_file, sep="\t", index=False)

        if d_thresh != 0:
            restricted_scaled_df = scaled_df[scaled_df.SCALED_DISTANCE > d_thresh]
            restricted_scaled_df.to_csv(full_out_file, sep="\t", index=False)
            plot_distances(restricted_scaled_df, out_prefix + "_restricted")

        # Get top 100 hits by distance.
        top100_df = scaled_df.head(100)
        top100_df.to_csv(top_out_file, sep="\t", index=False)

        # Plotting - only plots top 30k hits as browsers can't handle more.
        plot_distances(scaled_df, out_prefix + "_full")
        plot_distances(top100_df, out_prefix + "_top")

    print("Complete at: " + timeString() + ".")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_pre", required=True)
    parser.add_argument("-d", "--dthresh", dest="d_thresh", required=False, default=0, type=float)
    args = parser.parse_args()

    main(args.inp_file, args.out_pre, args.d_thresh)
