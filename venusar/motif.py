#!/usr/bin/env python3
"""
motif.py contains motif object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality
"""

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

    # Open provided motif file
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

                if (base_element.valid_flag):
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
        note: file format, can technically separate by separated by arbitrary
        strings of any whitespace characters
        (space, tab, newline, return, formfeed). Only space tested.

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
                # remove commas, brackets, and whitespace on far left and right
                line = line.strip().replace('[', '').replace(']', '').replace(',', '')
                if line != "":
                    line = line.split()
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