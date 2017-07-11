#!/usr/bin/env python3
"""
For a VAMPIRE'd VCF file, calculate a distance score for each motif delta score,
(optional) loci score, and gene expression score set for each variant. These values
will be plotted in 2 or 3 dimensional space. Additionally, simple bed-like files
are provided as output. One contains the scores for all motif, loci score, and
gene expression score sets for all variants. An optional second contains only
those sets that meet the distance score threshold as defined by the user.
A third will report the sets for the top 100 scaled-distance scores.

Usage: summarize.py -i <input.vcf> -o <output> [OPTIONS]

Args:
    -i (str): Path to sorted variant file to process.
    -o (str): Prefix for output files.
    -d (float, optional): Scaled distance magnitude threshold that must be met for variants/genes to be reported to
        output. Default is 0, so all set distances will be reported. A magnitude of 0.7 usually yields results pretty
        different from the norm for 3 dimensional data.
    -m (str): Path to motif file with counts.
    -a (float, optional) <0.25>: Background probability for A nucleotides. If not provided then all are assumed to
        be equally likely (all are 0.25).
    -c (float, optional) <0.25>: Background probability for C nucleotides.
    -g (float, optional) <0.25>: Background probability for G nucleotides.
    -t (float, optional) <0.25>: Background probability for T nucleotides.
    -pc (float, optional) <0.1>: Pseudocounts value to be added to all positions of the motif frequency matrix before
        calculating the probability matrix.
    -p (int, optional) <1>: Processor cores to utilize. Will decrease computation time linearly.
"""
import argparse
import sys
import time
from math import sqrt
from math import log10
from multiprocessing.dummy import Pool as ThreadPool
import itertools
from timeit import default_timer as timer

from Bio import motifs
import pandas as pd
import plotly
import plotly.graph_objs as go
import numpy as np
from pytfmpval import tfmp

import utils

# YYY-JA 04/24/2017 - Hate making yet another variant (Ha) of this class.
# Move to utils.py and make uniform across modules.


class Variant(object):
    """
    Use to process and handle variant records from a VCF more easily. Create from line of VCF file.
    """

    def __init__(self, line):
        self.line_list = line.strip().split("\t")
        self.pos = utils.Position(self.line_list[0], int(self.line_list[1]),
                                  (int(self.line_list[1]) + len(self.line_list[3])))
        self.ref_allele = self.line_list[3]
        self.var_allele = self.line_list[4]
        self.iden = self.line_list[2]
        self.orig_line = line.strip()
        self.info_fields = self.line_list[7].split(";")
        (self.common_var_samples, self.commons_samples, self.motif_fields, self.exp_fields, self.act_fields,
            self.genes) = self.parse_info_fields()
        self.motif_scores = self.get_motif_scores()  # Get var, ref, and diff motif scores and place into a dict.
        # Parse all combos of gene expression and loci data, determine if peaks are binary or quantitative.
        self.var_sample_data, self.binary = self.get_var_sample_data()

        self.output = self.get_variant_summary()  # Get output lines as a list of lists.

        if self.common_var_samples is not None:  # Should never evaluate to False.
            self.num_com_var_samps = len(self.common_var_samples)
        else:
            self.num_com_var_samps = 0

        if self.common_samples is not None:  # Should never evaluate to False.
            self.num_com_samps = len(self.common_samples)
        else:
            self.num_com_samps = 0

        self.var_recurrence_dec = self.num_com_var_samps / self.num_com_samps
        self.var_recurrence_frac = str(self.num_com_var_samps) + "/" + str(self.num_com_samps)

    def parse_info_fields(self):
        """
        Get names of samples containing variant and motif INFO fields from a
        variant record's INFO fields.

        Args:
            self (Variant): Variant object.

        Returns:
            common_samples (dict of tuple): List of samples that had variant
                called and have loci and expression data.
            motif_fields (list of str): List of INFO fields for variant that
                contain MOTIF related information.
            exp_fields (list of str): List of INFO fields for variant that
                contain Expression related information.
            act_fields (list of str): List of INFO fields for variant that
                contain Activity/Loci related information.
            genes (list of str): List of INFO fields for variant that
                contain MOTIF related information.
        """
        act_var_samples = None
        exp_var_samples = None
        act_samples = None
        exp_samples = None

        common_var_samples = []
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
                    exp_var_samples = data.split(",")
                if name == "EXPV" or name == "EXPR":
                    exp_samples = data.split(",")
            elif name.startswith("GENE"):
                genes = data.split(',')
            elif name.startswith("LOCI") or name.startswith("SAMP"):
                act_fields.append(field)
                # Get variant samples with locus activity data.
                if name == "SAMPSV":
                    act_var_samples = data.split(",")
                if name == "SAMPSV" or name == "SAMPSR":
                    act_samples = data.split(",")

        if act_var_samples:
            common_var_samples = compare_samples(exp_var_samples, act_var_samples)
        else:
            common_var_samples = exp_var_samples

        if act_samples:
            common_samples = compare_samples(act_samples, exp_samples)
        else:
            common_samples = exp_samples

        return (common_var_samples, common_samples, motif_fields, exp_fields, act_fields, genes)

    def get_motif_scores(self):
        """
        Returns the difference between reference and variant scores for each
        motif the variant matches as a dictionary of:
            {motif_name: (variant - reference log likelihood ratios)}.

        Returns:
            motifs (dict of floats): {motif_name: (variant_score, ref_score,
                               (variant - reference log likelihood ratios))}
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

        try:
            for i in names:
                idx = names.index(i)
                # Store scores for variant.
                motifs[i] = [float(var_scores[idx]), float(ref_scores[idx]),
                             float(var_scores[idx]) - float(ref_scores[idx])]
        except NameError:
            print("Motif names not found in file. Motif analysis must be performed first. Exiting.")
            sys.exit()

        return motifs

    def get_var_sample_data(self):
        """
        Parses and returns the gene expression and loci activity scores for variant samples.

        Returns:
            sample_data (list of lists of str):
                [[sample, loci1, sample act_z_score, gene1, sample exp_z_score],
                [sample, loci1, sample act_z_score, gene2, sample exp_z_score]...]
            binary (boolean): True if loci peaks are binary, False if quantitative.
        """

        samples = self.common_var_samples
        act_data = [x.split("=") for x in self.act_fields]
        exp_data = [x.split("=") for x in self.exp_fields]
        genes = self.genes

        sample_data = []
        binary = False

        if act_data:  # Handle if no activity data is used.
            for x in act_data:
                # Get loci id.
                if x[0] == "LOCIID":
                    loci = x[1].split(",")
                # Get each set of z-scores for each loci.
                elif x[0] == "LOCIVZ":
                    act_z_scores = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]
                # Get all variant samples that have the given peak if binary.
                elif x[0] == "LOCIVPK":
                    binary = True
                    peak_samples = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]

            # Create dict for loci - {loci_id: [peak_samples]}
            if binary:
                loci_data = {k: v for k, v in zip(loci, peak_samples)}
            # Or - {loci_id: [act_z_scores]}
            else:
                loci_data = {k: v for k, v in zip(loci, act_z_scores)}

        for x in exp_data:
            # Get each set of z-scores or fold-changes for each gene.
            if x[0] == "EXPVZ":
                gene_scores = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]
            elif x[0] == "EXPVFC":
                gene_scores = [x.strip("(").strip(")").split(',') for x in x[1].split("),(")]

        # Create dict for genes - {gene: [gene_scores]}
        gene_data = {k: v for k, v in zip(genes, gene_scores)}

        # Create list of lists containing all sets of relevant data combinations for all samples with a variant.
        # [[sample, loci1, sample act_data, gene1, sample exp_score],
        #  [sample, loci1, sample act_data, gene2, sample exp_score]...]
        if act_data:
            for s in samples:  # First iterate through all samples and get their expression and activity indices.
                e_idx = samples[s][0]
                a_idx = samples[s][1]

                for l in loci_data:  # Iterate through each locus for given sample.
                    samp_act_data = loci_data[l][a_idx]

                    for g in gene_data:  # And now each gene for each locus.
                        samp_exp_data = gene_data[g][e_idx]

                        sample_data.append([s, l, samp_act_data, g, samp_exp_data])
        else:  # Only return expression info if no activity data.
            for s in samples:  # First iterate through all samples and get their expression and activity indices.
                e_idx = samples[s][0]

                for g in gene_data:  # Get each gene for sample.
                    samp_exp_data = gene_data[g][e_idx]

                    sample_data.append([s, g, samp_exp_data])

        return (sample_data, binary)

    def get_variant_summary(self):
        """
        Create a list of summary output fields for a given Variant as well as
        its distance metrics for plotting.

        Returns:
            output_fields (list of lists of str): List of lists of output
                fields for the Variant.
        """
        var_info = [self.pos.chrom, self.pos.start, self.ref_allele, self.var_allele]
        motif_info = self.motif_scores
        var_sample_info = self.var_sample_data

        output_fields = []

        for m in motif_info:
            m_score = motif_info[m]  # Motif delta score for var vs ref.
            var_score, ref_score, diff_score = (m_score[0], m_score[1], m_score[2])
            for s in var_sample_info:
                if len(s) == 5:
                    sample = s[0]
                    loci = s[1]
                    loci_data = s[2]
                    gene = s[3]
                    exp_data = s[4]

                    if self.binary is False:
                        # To be used for plotting later.
                        dist_metrics = [(float(diff_score)), float(loci_data), float(exp_data)]
                        dist_score = calc_3d_dist(dist_metrics)
                        # Create complete list of output fields.
                        output_fields.append(var_info + [m] + [var_score] + [ref_score] + [diff_score] + [sample] +
                                             [loci] + [loci_data] + [gene] + [exp_data] + [dist_score] +
                                             [self.var_recurrence_frac] + [self.var_recurrence_dec])
                    else:  # TODO - How to determine "UNIQUE" and "UNIQUELY ABSENT" peaks?
                        dist_metrics = [(float(diff_score)), float(exp_data)]  # To be used for plotting later.
                        dist_score = calc_2d_dist(dist_metrics)
                        # Create complete list of output fields.
                        output_fields.append(var_info + [m] + [var_score] + [ref_score] + [diff_score] + [sample] +
                                             [loci] + [get_peak_status(loci_data)] + [gene] + [exp_data] +
                                             [dist_score] + [self.var_recurrence_frac] + [self.var_recurrence_dec])
                if len(s) == 3:  # If no loci activity data present.
                    sample = s[0]
                    gene = s[1]
                    exp_data = s[2]
                    dist_metrics = [(float(diff_score)), float(s[2])]  # To be used for plotting later.
                    dist_score = calc_2d_dist(dist_metrics)
                    # Create complete list of output fields.
                    output_fields.append(var_info + [m] + [var_score] + [ref_score] + [diff_score] + [sample] +
                                         [gene] + [exp_data] + [dist_score] + [self.var_recurrence_frac] +
                                         [self.var_recurrence_dec])

        return output_fields


def compare_samples(exp_samples, act_samples):
    """
    Compare two lists of samples and return common samples as a dict of format
    {sample_name: (list1_index, list2_index)} so data may be found appropriately.

    Args:
        exp_samples (list of str): List of variant samples with expression data.
        act_samples (list of str): List of variant samples with activity data.

    Return:
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


def get_peak_status(peak):
    """
    Return "Present" for 1 or "Absent" for 0.
    """

    if peak == 0:
        result = "ABSENT"
    elif peak == 1:
        result = "PRESENT"
    else:
        print("Binary peaks suspected, but 0 or 1 not present. Returning 'NA'.")
        result = "NA"

    return result


def calc_3d_dist(score_array):
    """
    Return a distance for a motif log-likelihood ratios, activity score, and gene expression score.
    """

    return sqrt((score_array[0] ** 2) + (score_array[1] ** 2) + (score_array[2] ** 2))


def calc_2d_dist(score_array):
    """
    Return a distance for a motif log-odds ratios and gene expression z-score.
    """

    return sqrt((score_array[0] ** 2) + (score_array[1] ** 2))


def check_peak_pres_binary(var_line):
    """
    Check if loci data is present and if it's binary or quantitative.

    Returns:
        present (boolean): True if loci peak data is present.
        binary (boolean): True if binary (0/1) loci data. False if quantitative.
    """
    binary = False
    present = False
    if len(var_line) == 11:
        loci_data = var_line[7]
        present = True

    if present:
        if all(type(x) is int for x in loci_data):
            binary = True

    return (present, binary)


def get_df_header(variant):
    """
    Use a Variant object to determine what the column headings of the dataframe should be. Adjusts header based on
    info fields found in variant record.

    Args:
        variant (Variant): Variant object to use to determine df column headers.

    Returns:
        header (list of string): List of column headers to be used by dataframe.
    """
    # These will always be present.
    header = ['CHR', 'POS', 'REF', 'ALT', 'MOTIF', 'VAR_SCORE', 'REF_SCORE', 'VAR-REF_SCORE', 'SAMPLE']

    act_fields = variant.act_fields
    exp_fields = variant.exp_fields

    if variant.act_fields:
        header.append("LOCIID")
        for field in act_fields:
            name = field.split('=')[0]
            if name == "LOCIVZ":
                header.append("LOCI_ZSCORE")
                break
            elif name == "LOCIVPK":
                header.append("PEAK")
                break
            elif name == "LOCIVFC":
                header.append("LOCI_LOG2_FC")
                break
    else:
        header.append("GENE")
        for field in exp_fields:
            name = field.split('=')[0]
            if name == "EXPVZ":
                header.append("LOCI_ZSCORE")
                break
            elif name == "LOCIVPK":
                header.append("PEAK")
                break
            elif name == "LOCIVFC":
                header.append("LOCI_LOG2_FC")
                break

    return header


def scale_and_frame(all_output, present, binary):
    """
    Return dataframe after adding scaled distance score to each row. Scaling done by dividing with max value of each
    metric, summation, and taking square root. Make things much easier to plot and handle. If loci data is presented
    and/or binary, will adjust the dataframe as appropriate.

    Args:
        all_output (list of lists): Each list contains a set of motif, expression, activity data.
            [chr, pos, ref, alt, motif, var score, ref score, diff score, sample, loci1, sample act_score,
             gene1, sample exp_score, distance]

    Return:
        df (pandas Dataframe): Each row contains a set of motif, expression, activity data.
            [chr, pos, ref, alt, motif, var score, ref score, diff score, sample, loci1, sample act_score,
             gene1, sample exp_score, distance, scaled_distance]
    """

    if present and not binary:
        df = pd.DataFrame(all_output, columns=['CHR', 'POS', 'REF', 'ALT', 'MOTIF', 'VAR_SCORE', 'REF_SCORE',
                                               'VAR-REF_SCORE', 'SAMPLE', 'LOCIID', 'ACT_ZSCORE', 'GENE',
                                               'EXP_ZSCORE', 'DISTANCE', 'RECURRENCE', '%RECURRENCE'])
        df[['VAR-REF_SCORE', 'VAR_SCORE', 'REF_SCORE',
            'ACT_ZSCORE', 'EXP_ZSCORE', 'DISTANCE']] = df[['VAR-REF_SCORE', 'VAR_SCORE', 'REF_SCORE', 'ACT_ZSCORE',
                                                           'EXP_ZSCORE', 'DISTANCE']].apply(pd.to_numeric)

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

    if present and binary:
        df = pd.DataFrame(all_output, columns=['CHR', 'POS', 'REF', 'ALT', 'MOTIF', 'VAR_SCORE', 'REF_SCORE',
                                               'VAR-REF_SCORE', 'SAMPLE', 'LOCIID', 'PEAK', 'GENE',
                                               'EXP_ZSCORE', 'DISTANCE'])
        df[['VAR-REF_SCORE', 'VAR_SCORE', 'REF_SCORE',
            'EXP_ZSCORE', 'DISTANCE']] = df[['VAR-REF_SCORE', 'VAR_SCORE', 'REF_SCORE',
                                             'EXP_ZSCORE', 'DISTANCE']].apply(pd.to_numeric)

        # Get scaling factors.
        motif_score_scale = max(df['VAR-REF_SCORE'] ** 2)
        gene_score_scale = max(df.EXP_ZSCORE ** 2)

        # Scale the distance and create new column.
        df['SCALED_DISTANCE'] = np.sqrt(((df['VAR-REF_SCORE'] ** 2) / motif_score_scale) +
                                        ((df.EXP_ZSCORE ** 2) / gene_score_scale))

        # Sort by distance.
        df.sort_values('SCALED_DISTANCE', ascending=False, inplace=True)

    return df


def get_pwms(motifs_file, background, pc=0.1):
    """
    Get PWMs for each motif in the motif_file.

    Args:
        motifs_file (str):
            Path to file containing motifs in JASPAR format as count matrices.
        pc (float):
            Pseudocount value to add to each position.
        background (list of float):
            A, C, T, and G background frequencies.

    Return:
        pwms (dict):
            {motif_name: pwm R matrix}
    """
    pwms = {}

    # Parse motifs using Biopython.
    fh = open(motifs_file)
    for m in motifs.parse(fh, "jaspar"):
        pfm = m.counts.normalize(pseudocounts=pc)    # Create frequency matrix.
        pwm = pfm.log_odds(background)              # Calculate to log likelihoods vs background.
        # Create matrix string from motif pwm.
        mat = pwm[0] + pwm[1] + pwm[2] + pwm[3]
        mat = [str(x) for x in mat]
        mat = " ".join(mat)
        pwms[m.name] = mat
    fh.close()

    return pwms


def find_pval(matrices, m_score, r_bg):
    """
    Utilize pytfmpval to calculate the p-values for a score and motif given
    the nucleotide background frequencies and PWM list.

    Args:
        matrices (dict of R matrix):
            Dict of PWMs corresponding to each motif. {motif name: R PWM}
        m_score (tuple):
            Tuple of (score, motif).
        r_bg (R FloatVector):
            Vector containing background nucleotide frequencies [A, C, G, T]
    """
    score = m_score[0]
    motif = m_score[1]
    matrix = matrices[motif]

    start = timer()
    mat = tfmp.read_matrix(matrix, r_bg)  # Create the pytfmpval Matrix object.
    pval = tfmp.score2pval(mat, score)  # Actually calculate p-value from score.

    if pval[0] == 0:
        pval[0] = 0.000000000000001  # So log can still be taken.
    end = timer()

    print(motif + ": " + str(end - start))

    return float(pval[0])


def get_pvals(df, bg, matrices, proc):
    """
    Utilize the pytfmpval package to calculate p-values for both the variant and reference scores for a given
    variant and motif. Add them to the DataFrame as additional columns. Also add the log2 fold-change between them.

    Args:
        df (pandas DataFrame):
            Dataframe containing the motif scores for the variant and reference samples.
        bg (list of float):
            A, C, T, and G background frequencies.
        matrices (dict of R matrices):
            {motif name: R pwm matrix}
        proc (int):
            Processor cores to be used.

    Return:
        df (pandas DataFrame):
            Updated DataFrame with p-values.
    """

    ref_scores = [x for x in zip(df.REF_SCORE, df.MOTIF)]
    var_scores = [x for x in zip(df.VAR_SCORE, df.MOTIF)]

    with ThreadPool(proc) as pool:
        var_pvals = pool.starmap(find_pval, zip(itertools.repeat(matrices), var_scores,
                                                itertools.repeat(bg)))
    with ThreadPool(proc) as pool:
        ref_pvals = pool.starmap(find_pval, zip(itertools.repeat(matrices), ref_scores,
                                                itertools.repeat(bg)))

    # Insert columns for p-values.
    df.insert(6, 'VAR_PVAL', var_pvals)
    df.insert(8, 'REF_PVAL', ref_pvals)
    df[['VAR_PVAL', 'REF_PVAL']] = df[['VAR_PVAL', 'REF_PVAL']].apply(pd.to_numeric)

    return df


def plot_3d_dist(df, out_prefix):
    """
    Plot distance scores calculated for the delta var/ref motif log-odds ratios,
    activity z-score, and gene expression z-score sets for each variant.

    Args:
        df (pandas Dataframe) = Dataframe containing all variant, distance
            metric, and sample info. One record per row.
        out_prefix (str) = Prefix to use for plot outputs.
    """

    # Take top 30k hits only, plotly can't really handle more for 3D plots.
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
            colorscale=[[0, 'rgb(188, 191, 196)'], [0.5, 'rgb(188, 191, 196)'], [1, 'rgb(170, 1, 15)']],   # colorscale
            colorbar=dict(
                title='Scaled Distance'
            ),
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
                title='Var/Ref Motif Log-Likelihood Ratio Difference',
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
    plotly.offline.plot(
        fig, filename=out_prefix + '.html',
        auto_open=False, image="png", image_filename=out_prefix
    )


def plot_2d_vol(df, out_prefix):
    """
    Create volcano plot for samples compared to a control group.

    X-axis is gene expression fold-change as compared to control group. Variants in the
    control group are obviously ignored and excluded from output.

    Args:
        df (pandas Dataframe) = Dataframe containing all variant, distance
            metric, and sample info. One record per row.
        out_prefix (str) = Prefix to use for plot outputs.
    """

    # Take top 30k hits only, plotly can't really handle more for 3D plots.
    if len(df) > 30000:
        df = df.head(30000)

    info = list(zip(df.SAMPLE, df.MOTIF, df.GENE))
    info_list = ["Sample: " + x[0] + ", Motif: " + x[1] + ", Gene: " + x[2] for x in info]

    trace1 = go.Scatter(
        name="Distances",
        x=df['VAR-REF_SCORE'],
        y=-log10(df.VAR_PVAL),
        z=df.EXP_ZSCORE,
        hoverinfo="x+y+text",
        text=info_list,
        mode='markers',
        marker=dict(
            size=4,                 # Using gene expression as color scale values for now.
            color=df.GROUP,                # set color to an array/list of desired values
            colorscale=[[0, 'rgb(188, 191, 196)'], [0.5, 'rgb(188, 191, 196)'], [1, 'rgb(170, 1, 15)']],   # colorscale
            colorbar=dict(
                title='Scaled Distance'
            ),
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
        title='Volcano Plot for Individual Variant Events',
        scene=dict(
            xaxis=dict(
                title='Variant Motif p-value',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            ),
            yaxis=dict(
                title='Gene Expression Fold-Change',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            ),
        )
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(
        fig, filename=out_prefix + '.html',
        auto_open=False, image="png", image_filename=out_prefix
    )


def plot_2d_cor(df, out_prefix):
    """
    Create volcano plot for samples compared to a control group.

    X-axis is gene expression fold-change as compared to control group. Variants in the
    control group are obviously ignored and excluded from output.

    Args:
        df (pandas Dataframe) = Dataframe containing all variant, distance
            metric, and sample info. One record per row.
        out_prefix (str) = Prefix to use for plot outputs.
    """

    # Take top 30k hits only, plotly can't really handle more for 3D plots.
    if len(df) > 30000:
        df = df.head(30000)

    info = list(zip(df.SAMPLE, df.MOTIF, df.GENE))
    info_list = ["Sample: " + x[0] + ", Motif: " + x[1] + ", Gene: " + x[2] for x in info]

    trace1 = go.Scattergl(
        name="Distances",
        x=df.EXP_LOG2_FC,
        y=-log10(df.VAR_PVAL),
        z=df.EXP_ZSCORE,
        hoverinfo="x+y+text",
        text=info_list,
        mode='markers',
        marker=dict(
            size=4,                 # Using gene expression as color scale values for now.
            color=df.SCALED_DISTANCE,                # set color to an array/list of desired values
            colorscale=[[0, 'rgb(188, 191, 196)'], [0.5, 'rgb(188, 191, 196)'], [1, 'rgb(170, 1, 15)']],   # colorscale
            colorbar=dict(
                title='Scaled Distance'
            ),
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
        title='Volcano Plot for Individual Variant Events',
        scene=dict(
            xaxis=dict(
                title='Var/Ref Motif Log-Likelihood Ratio Difference',
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
    plotly.offline.plot(
        fig, filename=out_prefix + '.html',
        auto_open=False, image="png", image_filename=out_prefix
    )


def main(vcf_file, out_prefix, d_thresh, motif_file, background, pc, proc):
    """
    Summarize each variant/motif/gene/sample/enhancer set, calculating p-values, distance scores, and plots.

    Args:
        vcf_file (str):
            Path to sorted variant file to process.
        out_prefix (str):
            Prefix to be used for output files.
        d_thresh (float):
            Distance threshold to be used for more restricted plotting and reporting.
        motif_file (str, optional):
            Path to motif file.
        background (list of float):
            Background frequency of A, C, G, and T. Must sum to 1.
        pc (float):
            Pseudocount value to be used for PWM generation.
        proc (int):
            Number of processor cores to utilize.
    """

    with open(vcf_file) as f:

        line = f.readline().strip()
        now = time.strftime("%c")

        command = ('##venusar=<ID=summary,Date="' + now + '",CommandLineOptions="--input ' + vcf_file +
                   ' --output ' + out_prefix + ' --dthresh ' + str(d_thresh) + ' --motif ' + motif_file + '">')
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

        # Check if loci data is present and if it's binary.
        present, binary = check_peak_pres_binary(all_output[0])

        # Scale distance metrics.
        print("Calculating distance metrics.")
        scaled_df = scale_and_frame(all_output)

        # Get p-values if motif file provided.
        if motif_file:
            print("Calculating p-values. This may take some time.")
            bio_background = {'A': background[0], 'C': background[1], 'G': background[2], 'T': background[3]}
            pwms = get_pwms(motif_file, bio_background, pc)
            scaled_df = get_pvals(scaled_df, background, pwms, proc)

        print("Creating output files and plots.")
        scaled_df.to_csv(full_out_file, sep="\t", index=False)

        if d_thresh != 0:
            restricted_scaled_df = scaled_df[scaled_df.SCALED_DISTANCE > d_thresh]
            restricted_scaled_df.to_csv(rest_out_file, sep="\t", index=False)
            plot_3d_dist(restricted_scaled_df, out_prefix + "_restricted")

        # Get top 100 hits by distance.
        top100_df = scaled_df.head(100)
        top100_df.to_csv(top_out_file, sep="\t", index=False)

        # Plotting - only plots top 30k hits as browsers can't handle more.
        plot_3d_dist(scaled_df, out_prefix + "_full")
        plot_3d_dist(top100_df, out_prefix + "_top")

    print("Complete at: " + utils.timeString() + ".")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", "--input", dest="inp_file", required=True)
    parser.add_argument("-o", "--output", dest="out_pre", required=True)
    parser.add_argument("-d", "--dthresh", dest="d_thresh", required=False, default=0, type=float)
    parser.add_argument("-m", "--motif", dest="motif_file", required=False, default=None)
    parser.add_argument("-a", "--a_freq", dest="a_freq", required=False, default=0.25, type=float)
    parser.add_argument("-c", "--c_freq", dest="c_freq", required=False, default=0.25, type=float)
    parser.add_argument("-g", "--g_freq", dest="g_freq", required=False, default=0.25, type=float)
    parser.add_argument("-t", "--t_freq", dest="t_freq", required=False, default=0.25, type=float)
    parser.add_argument("-pc", "--pseudocounts", dest="pseudocounts",
                        required=False, default=0.1, type=float)
    parser.add_argument("-p", "--processors", dest="processors", required=False, default=1, type=int)
    args = parser.parse_args()

    if sum([args.a_freq, args.c_freq, args.g_freq, args.t_freq]) == 1.0:
        bp = [args.a_freq, args.c_freq, args.g_freq, args.t_freq]
    else:
        print("Background frequencies must equal 1. Check input parameters, exiting.")
        sys.exit()

    main(args.inp_file, args.out_pre, args.d_thresh, args.motif_file, bp, args.pseudocounts, args.processors)
