#!/usr/bin/env python3
"""
motif.py contains motif object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality
"""

import sequence    # need sub_from_*, crop_from_*
                   # (maybe later need int func and split wing/core processing)
                   # then would want base sequence object
from math import log2

class motifArray:
    """this class is the array of motif Elements"""

    def __init__(self):
        self.max_positions = 0    # largest positions value for all Elements
        self.max_positions_index = -1    # index of maximum position length
        self.motifs = []
        return

    def add_motif(self, motif_element2add):
        """add passed element to the set"""
        self.motifs.append(motif_element2add)
        if (motif_element2add.positions > self.max_positions):
            self.max_positions = motif_element2add.positions
            self.max_positions_index = len(self.motifs) - 1

        return

    def delete_motif(self, motif_index):
        """remove passed element index from the set"""
        if (motif_index < 0 | motif_index > len(self.motifs)):
            return

        reset_max_length = False
        if len(self.motifs) > 1:
            if motif_index == self.max_positions_index:
                reset_max_length = True
        else:
            self.__init__()    # reset to defaults
            return

        # get to here then know started with multi-element set
        del self.motifs[motif_index]

        if reset_max_length:
            self.search_set_max_position()

        return

    def join(self, motif_array2):
        """join another motifArray object to this motifArray object"""

        # check validity of the items to be joined
        if type(motif_array2) is not motifArray:
            return
        if motifArray.length() == 0:
            return

        # empty current, non-empty new --> replace
        if self.length() == 0:
            self = motif_array2    # allowed in python?
            return

        # comparison elements
        current_max = [self.max_positions, self.max_positions_index]
        joined_max = [motif_array2.max_positions, motif_array2.max_positions_index]
        original_length = self.length()

        # -- compare and join
        self.motifs.append(motif_array2.motifs)
        if (joined_max[0] > current_max[0]):
            self.max_positions = joined_max[0]
            self.max_positions_index = original_length + joined_max[1]

        return

    def length(self):
        """return length of the set of motifs for the array"""
        return (len(self.motifs))

    def motif_match(self, baseline_p, ref_seq, var_seq, wing_l):
        """
        Takes a reference and variant sequence string, then checks for matches
        in the motif set. Outputs match score for both ref seq and var seq
        for any case where either matches above the motif element threshold.
        (former name: match_motifs, wanted grouped in class function list)

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            ref_seq = Reference sequence represented as a string (ACGT only)
            var_seq = Variant sequence represented as a string (ACGT only)
            wing_l = Integer length of sequence of bases flanking the variant
                (generally >= self.max_positions, should match value used to
                create ref_seq and var_seq)

        Returns: A list of Motif_match objects
        """

        # scoring returns: motif match array with the form
        #    (id, name, threshold, max_score, match_seq)
        scored = self.motif_scores(baseline_p, var_seq, wing_l)
        r_scored = self.motif_scores(baseline_p, ref_seq, wing_l)

        matches = []

        # generate list of motifs that matched var seq to compare ref seq matches
        for idx in range(len(scored)):
            (iden, name, th, max_score, match_seq) = scored[idx]
            (rid, rname, rth, rmax_score, rmatch_seq) = r_scored[idx]
            if (max_score >= th or rmax_score >= rth):
                # XXX:QQQ: there needs to be a margin of error for th != rth; using 1%. appropriate?
                if (iden != rid or name != rname or abs((th - rth) / th) > .01):
                    print(("***ERROR*** matching motifs to varseq and refseq desynced\n" +
                          iden + " != " + rid +
                          " or " + name + " != " + rname +
                          " or " + th + " != " + rth))
                # tup = (id, name, max_score, match_seq, rmax_score, rmatch_seq)
                match = Motif_match(name, max_score, rmax_score)
                matches.append(match)

        # output variant seq matches vs ref seq matches
        """for m in matches:
            line = "Variant matched motif "+m[1]+" ("+m[0]+") with score "
            line+= str(m[2])[:6]+" at seq: "+m[3]+"\nReference matched motif "+m[1]
            line+= "           with score "+str(m[4])[:6]+" at seq: "+m[5]
            #debug line+= "\n\tRefseq: "+ref_seq+"\n\tVarseq: "+var_seq
            print(line,file=output_f)"""    # WARNING: debug using undefined global

        return matches

    def motif_scores(self, baseline_p, sequence_str, wing_l, normalize):
        """
        Calculate if any motifs in the motif list match the given sequence.
        Requires that no motif have a length of 0.
        (former name: score_motifs, wanted grouped in class function list)

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')
            wing_l = length of sequence of bases flanking the variant
            normalize = boolean, if true then divide score by motif length
                concern: if not normalized, longer motifs can generate higher
                scores by length not true matches.

        Returns:
            matches = array of scores of the form:
                ( id,  name,  threshold,  max_score,  match_seq,  motif_index )

        """
        # list of matches between motifs and sequence
        scores = []

        for motif_index in range(len(self.motifs)):
            #### (iden, name, thresh, p_matrix) = tup
            motif_element = self.motifs[motif_index]

            # -- prepare the motif element and sequence string for scoring
            if not motif_element.valid_flag:
                continue

            if motif_element.matrix_type == 0:
                # need to calculate base occurence frequency by position. why
                # didn't it do this during insertion? will call many times here
                print("motif_scores: probabilities not previously calculated for " +
                    motif_element.print_str + ". Likely slower to calculate here. Fix")
                motif_element.calculate_probabilities()

            # trim flanking bases to one less than motif length
            # number of bases to trim
            trim_amount = wing_l - motif_element.positions + 1
            trim_seq = sequence.crop_from_left(sequence_str, trim_amount)
            trim_seq = sequence.crop_from_right(trim_seq, trim_amount)

            # -- actually do the scoring
            # highest match score for this motif (based on starting position)
            max_score = float("-inf")
            # sequence that matched motif best
            match_seq = ""

            # check match to all positions where motif overlaps with variant
            for pos in range(motif_element.positions):
                # check match starting at position
                # (score_motif will stop after the length of the motif)
                pos_score = motif_element.motif_score(baseline_p, trim_seq[pos:])
                if pos_score > max_score:
                    max_score = pos_score
                    match_seq = trim_seq[pos:pos + motif_element.positions]

            if normalize:
                max_score = max_score / motif_element.positions

            # debug print("Max match score:"+str(max_score)+" for motif "+name+" and
            # sequence "+match_seq+".")
            tupl = (motif_element.id, motif_element.name,
                    motif_element.threshold, max_score, match_seq, motif_index)
            scores.append(tupl)

        return scores

    def process_local_env(self, baseline_p, matches, seq_element, co_binders_dict, v_seq, r_seq):
        """
        Finds GC content, weak homotypic matches, and motif matches for co-binding
            transcription factors. Co-binding tf currently not implemented.

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            matches = list of Motif_match objects, generated by motif_match
            seq_element = sequenceElement object v_seq and r_seq belong to
            co_binders = dictionary that lists co-binding transcription factors for
                a given transcription factor name [feature not enabled. Use None]
            v_seq = variant sequence array of characters (strings 'A','C','G', or 'T')
            r_seq = reference sequence array of characters (strings 'A','C','G', or 'T')

        Returns: Updates matches
        XXX: this function maybe isn't properly updating class elements;
            review call usage and code.
        XXX: WARNING: ERROR: this function is confusing class and other elements
            of the match
        """

        # Get GC content [QQQ: faster to convert to upper then count?]
        gc = v_seq.count('G') + v_seq.count('g')
        gc += v_seq.count('C') + v_seq.count('c')
        var_gc = gc / len(v_seq)

        gc = r_seq.count('G') + r_seq.count('g')
        gc += r_seq.count('C') + r_seq.count('c')
        ref_gc = gc / len(r_seq)

        # Get list of motifs that matched
        for matches_index in range(len(matches)):
            match = matches[matches_index]    # XXX
            match_index = match[5]
            if match_index < 0 or match[5] > self.length():
                print(("**Error** process_local_env unable to find motif" +
                    match[1] + ":id:" + match[0] ))
                continue;

            # motif info
            #self.motifs[match_index].matrix
            #self.motifs[match_index].threshold
            match_motif = self.motifs[match_index]

            # Find homotypic variant and reference matches
            vh_matches = []
            #rh_matches = []    # do not compute, same as vh_matches (XXX: test presumption)

            # Check each wing seperately,
            # homotypic matches should not overlap the variant
            # because variant and reference wings match only compute once
            # wing size = size of the motif matched?
            left_wing = sequence.sub_from_left(seq_element.seq_left_wing.seq, match_motif.positions)
            right_wing = sequence.sub_from_right(seq_element.seq_right_wing.seq, match_motif.positions)
            vh_matches = match_motif.ht_matches_in(baseline_p, left_wing)
            vh_matches += match_motif.ht_matches_in(baseline_p, right_wing)

            # Update the match information: XXX
            match.var_ht = vh_matches
            match.var_gc = var_gc
            #match.ref_ht = rh_matches    # do not compute because same as var
            match.ref_gc = ref_gc

        # Finding co-binders currently not implemented

        return matches    # XXX: wasn't returning anything at all before

    def search_set_max_position(self):
        """search elements in array to find and set max position information"""

        # set values to defaults
        self.max_positions = 0
        self.max_positions_index = -1

        # loop to set right value
        #for element in self.motifs:
        for index in range(0, (len(self.motifs) - 1)):
            if (self.motifs[index].positions > self.max_positions):
                self.max_positions = self.motifs[index].positions
                self.max_positions_index = index

        return


class motifElement:
    """
    This class creates a Motif object
    Each object defines Motif characteristics principally the
    position matrix, matrix size, and modifications to the matrix
    """

    def __init__(self):
        self.name = "undefined"
        self.id = "undefined"
        self.matrix = [[] for x in range(4)]    # placeholder
            # matrix columns are positions
            # matrix rows are for bases, order ACGT
        self.matrix_type = 0    # 0 = count (Default), 1 = probability
                                # set by calculate_probabilities()
        self.positions = 0      # number of positions for motif matrix
                                # set by check_valid()
        self.valid_flag = False  # set by check_valid()
            # QQQ: store both types of matrix? memory hit is likely small
        self.threshold = 0
        self.base_probability = [] * 4   # base pair probability
        self.pseudocount = 0    # count number to add before probability calc

    def calculate_probabilities(self):
        """
        Given a count motif matrix, return motif probabilities
            matrix columns are positions
            matrix rows are for bases, order ACGT    XXX|QQQ make sure notes match

        Each input matrix element represents the number of times that base
        appears at that position in that motif for base defined by row and
        position defined by column.
        Each occurrences count is divided by the total number of occurrences
        for that position to give the frequency of each base by position.
        The position columns sum to 1.

        Pseudocounts are added in to avoid having a zero value in matrices with
        small sample sizes; see motifElement.pseudocount
        """

        # only process count matrices that are valid
        if (self.matrix_type != 0 or not(self.valid_flag)):
            return

        # for clarity (also note: self.positions = length any matrix row)
        a = self.matrix[0]
        c = self.matrix[1]
        g = self.matrix[2]
        t = self.matrix[3]

        # QQQ: change to use numpy methods --> faster?
        # initialize output matrix of same size as input matrix
        new_m = [[0 for y in range(self.positions)] for x in range(4)]

        # divide each position by the sum of its column
        for pos in range(self.positions):
            total = float(a[pos]) + float(c[pos]) + float(g[pos]) + float(t[pos])
            total += 4 * self.pseudocount
            for base in range(4):
                new_m[base][pos] = (float(self.matrix[base][pos]) +
                     self.pseudocount) / total

        self.matrix = new_m
        self.matrix_type = 1

        return

    def check_valid(self):
        """
            sets valid_flag and positions for the motifElement
            if all matrix elements are equal size --> set valid_flag = True
            if any matrix row length different than base length --> = False
        """
        self.valid_flag = True
        for mIndex in range(len(self.matrix)):
            if (mIndex == 0):
                self.positions = len(self.matrix[mIndex])
            else:
                if (self.positions != len(self.matrix[mIndex])):
                    self.valid_flag = False
                    break

        return

    def clear(self):
        self = self.__init__()    # ZZZ: note: self = probably not necessary
        return

    def ht_matches_in(self, baseline_p, sequence_str):
        """
        Compute homotypic matches (scores above threshold) in the given sequence

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')
        Returns:

        see also: motifArray.motif_scores()
        """
        homotypic_matches = []
        for pos in range(len(sequence_str) - len(self.positions) + 1):
            pos_score = self.score_motif(baseline_p, sequence_str[pos:])
            # Return the homotypic match if it is above the threshold
            if pos_score > self.threshold:
                homotypic_matches.append(pos_score)
        return homotypic_matches

    def print_str(self):
        """return a string of the class object state"""
        print_string = self.name
        print_string += (":id:" + self.id)
        print_string += (":valid:" + self.valid_flag)
        print_string += (":matrix_type:" + self.matrix_type)
        print_string += (":positions:" + format(self.positions))
        print_string += (":threshold:" + format(self.threshold))
        print_string += (":pseudoCount:" + format(self.pseudocount))
        print_string += (":base_probability:" + format(self.base_probability))
        print_string += (":matrix:" + self.matrix)
        return print_string

    def score_base(self, probability, baseline_p):
        """
        compute the ... log likelihood ratio
        XXX: not in score motif says natural log but this is log2
        YYY: note this could be outside the class
        """
        if (probability == 0):
            # should mathematically be negative infinity
            # this should not occur with pseudocounts added in
            return -100
        else:
            return log2(probability / baseline_p)

    def score_motif(self, baseline_p, sequence_str):
        """
        Calculate match between probability matrix for motifElement and the
            given sequence_str starting at the beginning of the sequence
            and ending after the length of the motif.

        Args:
            p_matrix = probability matrix where rows are bases (A,C,G,T) and columns
                are probabilities (that sum to 1) at a given position. All rows should
                be the same length.
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')

        Returns:
            Match score between probability matrix and sequence_str
        """

        score = 0

        if len(sequence_str) < len(self.positions):
            print("Sequence shorter than probability matrix. Returning score of 0.")
            return score

        if not self.valid_flag or self.matrix_type != 1:
            return score

        for pos in range(len(self.positions)):
            # Match base
            if (sequence_str[pos] == 'A' or sequence_str[pos] == 'a'):
                # natural log of likelihood ratio    ### XXX but computes log2!
                score += self.score_base(self.matrix[0][pos], baseline_p[0])
            elif (sequence_str[pos] == 'C' or sequence_str[pos] == 'c'):
                score += self.score_base(self.matrix[1][pos], baseline_p[1])
            elif (sequence_str[pos] == 'G' or sequence_str[pos] == 'g'):
                score += self.score_base(self.matrix[2][pos], baseline_p[2])
            elif (sequence_str[pos] == 'T' or sequence_str[pos] == 't'):
                score += self.score_base(self.matrix[3][pos], baseline_p[3])
            # else is score = 0 therefore do nothing

        return score


class Motif_match:
    """
    This class could be rolled into motifElement
    Current theory; faster not to, like normalized database table set,
    not duplicating motifElement variables don't need and don't need to loop
    across the full motifArray set
    """
    def __init__(self, name, var_score, ref_score):
        # The transcription factor (tf) name
        self.name = name

        # tf motif match information (variant and reference sequence max log
        # likelihood score for the tf)
        self.var_score = var_score
        self.ref_score = ref_score

        # Does the motif match fall inside a chip peak for the same tf
        self.chip_match = False

        # Sequence environment data
        # variant and reference gc content
        self.var_gc = None
        self.ref_gc = None
        # variant and reference homotypic matches
        self.var_ht = []
        self.ref_ht = []    # XXX: drop see process_local_env notes
        # co-binding motifs (currently not implemented)
        self.cobinders = []

        self.motif_array_index = -1    # for corresponding motifArray object


# ______ START METHODS RELATED TO THE CLASS BUT NOT TIED TO IT _____
#    many could be made part of motifArray or motifElement

def get_motifs(motif_filename, pc, default_th, base_pr):
    """
    Read in and calculate probability matrices from a frequency matrix file.

    Args:
        motif_f = filename of the frequency matrix
        pc = pseudocount value (small, greater than zero) to be added to all
            positions in matrix
        default_th = default threshold value used if none is listed in the
            motif file (should be third part of header if present).
            If None is given, default threshold will be calculated for motifs.
        base_pr = probabilities as a float array of the form:
            [ PrA, PrC, PrG, PrT ]
            ** Currently unused by this method
                but could be used to return a pssm instead of a pwm

    Returns:
        MotifArray object = set of Motif sequences as motifElement objects
        Each object in the array defines Motif characteristics

        If reading multiple files use MotifArray join method on returned object


    XXX: QQQ: believe this should index motif length by TF.
        also modify wing length on individual not max for all
    """

    motif_set = motifArray()

    # Open and import motif file: note: with always closes open file
    with open(motif_filename) as file_handle:

        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n
        # arrays for C, G, T - each with same format as A
        # 1 blank line separates each individual TF position weight matrix
        # reading files assumes {header,A,C,G,T,blankline} order
        base_element = motifElement()

        # index for line you are looking at in the motif matrix set (not file)
        i = 0

        # Iterate through file line by line.
        for line in file_handle:

            # First line contains id and name
            if i == 0:
                # must first remove brackets and split on whitespace
                line = line.strip().strip('[').strip(']').split()

                base_element.id = line[0]
                base_element.name = line[1]
                    #print(format(i) + ":" + base_element.id + " " + base_element.name)    # DEBUG LINE
                # If threshold given, use it for this matrix. Else, use default.
                if (len(line) > 2):
                    base_element.threshold = float(line[2])
                else:
                    base_element.threshold = default_th

            # Order of position weight matrices is A,C,G,T. Remove brackets etc.
            elif i < 5:
                # must first remove brackets and split on whitespace
                #     technically should modify index by found ACGT
                    #print(i)    # DEBUG LINE
                line = line.strip().strip('A').strip('C').strip('G').strip('T').strip().strip('[').strip(']').split()
                base_element.matrix[i - 1] = line[0:]

            # add motif to list and continue (there are 2 newlines between motifs)
            else:
                # check validity and set values passed
                base_element.check_valid()
                base_element.base_probability = base_pr
                base_element.pseudocount = pc

                if base_element.valid_flag:
                    # calculate base occurence frequency by position
                    base_element.calculate_probabilities()

                # QQQ if invalid should it not include?
                # append the current TF matrix to the motif set
                motif_set.add_motif(base_element)

                # reset the tracker variables
                i = -1
                base_element.clear()
            i += 1

    return motif_set

def get_baseline_probs(baseline_f):
    """
    Read in baseline probabilities from a file name.

    Args:
        baseline_f a file containing a probability array of the form:
            [ PrA PrC PrG PrT ]
        Where PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)
        note1: the file can contain header lines marked by #
        note2: file format can technically be separated by arbitrary strings of
            any whitespace characters (space, tab, newline, return, formfeed).
            Only space tested during code development.
        note3: Separator must be whitespace, commas converted to empty string.
        note4: subsequent lines ignored once a properly formatted line is found

    Returns:
        Array with probabilities as a float array of the form:
        [ PrA, PrC, PrG, PrT ]
    """

    # Default baseline probability numbers (assumes all are equally likely)
    bp_array = [0.25, 0.25, 0.25, 0.25]

    # QQQ|YYY: note safer if check for invalid files everywhere
    with open(baseline_f) as f:
        try:
            for line in f:
                # remove header
                if line.startswith("#"):
                    continue
                # remove commas, brackets, and whitespace on far left and right
                line = line.strip().replace('[', '').replace(']', '').replace(',', '')
                if line != "":
                    line = line.split()
                    if (len(line) >= 4):
                        for idx in range(4):
                            bp_array[idx] = float(line[idx])
                        return bp_array
        except ValueError:
            print(("**ERROR** Baseline probability file incorrectly formatted.\n" +
                  "\tFile should contain only [ PrA PrC PrG PrT ] \n" +
                  "\tWhere PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)\n" +
                  "\tContinuing with default probabilities: " + format(bp_array)))
            return bp_array

        print(("**ERROR** Empty file found.\n" +
              "\tFile should contain only [ PrA PrC PrG PrT ] \n" +
              "\tWhere PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)\n" +
              "\tContinuing with default probabilities: " + format(bp_array)))
        return bp_array