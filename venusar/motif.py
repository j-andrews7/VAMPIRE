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
        self.max_positionsIndex = -1    # index of maximum position length
        self.motifs = []

    def addMotif(self, motifElementToAdd):
        """add passed element to the set"""
        self.motifs.append(motifElementToAdd)
        if (motifElementToAdd.positions > self.max_positions):
            self.max_positions = motifElementToAdd.positions
            self.max_positionsIndex = len(self.motifs) - 1

        return

    def deleteMotif(self, motifIndex):
        """remove passed element index from the set"""
        if (motifIndex < 0 | motifIndex > len(self.motifs)):
            return

        resetMaxLength = False
        if len(self.motifs) > 1:
            if motifIndex == self.max_positionsIndex:
                resetMaxLength = True
        else:
            self.__init__()    # reset to defaults
            return

        # get to here then know started with multi-element set
        del self.motifs[motifIndex]

        if resetMaxLength:
            self.searchSetMaxPosition()

        return

    def join(self, motifArray2):
        """join another motifArray object to this motifArray object"""

        # check validity of the items to be joined
        if type(motifArray2) is not motifArray:
            return
        if motifArray.length() == 0:
            return

        # empty current, non-empty new --> replace
        if self.length() == 0:
            self = motifArray2    # allowed in python?
            return

        # comparison elements
        currentMax = [self.max_positions, self.max_positionsIndex]
        joinedMax = [motifArray2.max_positions, motifArray2.max_positionsIndex]
        originalLength = self.length()

        # -- compare and join
        self.motifs.append(motifArray2.motifs)
        if (joinedMax[0] > currentMax[0]):
            self.max_positions = joinedMax[0]
            self.max_positionsIndex = originalLength + joinedMax[1]

        return

    def length(self):
        """return length of the set of motifs for the array"""
        return (len(self.motifs))

    def searchSetMaxPosition(self):
        """search elements in array to find and set max position information"""

        # set values to defaults
        self.max_positions = 0
        self.max_positionsIndex = -1

        # loop to set right value
        #for element in self.motifs:
        for index in range(0, (len(self.motifs) - 1)):
            if (self.motifs[index].positions > self.max_positions):
                self.max_positions = self.motifs[index].positions
                self.max_positionsIndex = index

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
        self.matrixType = 0     # 0 = count (Default), 1 = probability
                                # set by calculate_probabilities()
        self.positions = 0      # number of positions for motif matrix
                                # set by checkValid()
        self.validFlag = False    # set by checkValid()
            # QQQ: store both types of matrix? memory hit is likely small
        self.threshold = 0
        self.base_probability = [] * 4   # base pair probability
        self.pseudoCount = 0    # count number to add before probability calc

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
        small sample sizes; see motifElement.pseudoCount
        """

        # only process count matrices that are valid
        if (self.matrixType != 0 or not(self.validFlag)):
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
            total += 4 * self.pseudoCount
            for base in range(4):
                new_m[base][pos] = (float(self.matrix[base][pos]) +
                     self.pseudoCount) / total

        self.matrix = new_m
        self.matrixType = 1

        return

    def checkValid(self):
        """
            sets validFlag and positions for the motifElement
            if all matrix elements are equal size --> set validFlag = True
            if any matrix row length different than base length --> = False
        """
        self.validFlag = True
        for mIndex in range(len(self.matrix)):
            if (mIndex == 0):
                self.positions = len(self.matrix[mIndex])
            else:
                if (self.positions != len(self.matrix[mIndex])):
                    self.validFlag = False
                    break

        return

    def clear(self):
        self = self.__init__()    # ZZZ: note: self = probably not necessary

    def printStr(self):
        printStr = self.name
        printStr += (":id:" + self.id)
        printStr += (":valid:" + self.validFlag)
        printStr += (":matrixType:" + self.matrixType)
        printStr += (":positions:" + self.positions)
        printStr += (":threshold:" + self.threshold)
        printStr += (":pseudoCount:" + self.pseudoCount)
        printStr += (":base_probability:" + self.base_probability)
        printStr += (":matrix:" + self.matrix)


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

    motifSet = motifArray()

    # Open provided motif file
    with open(motif_filename) as fileHandle:

        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n
        # arrays for C, G, T - each with same format as A
        # 1 blank line separates each individual TF position weight matrix
        # reading files assumes {header,A,C,G,T,blankline} order
        baseElement = motifElement()

        # index for line you are looking at in the motif matrix set (not file)
        i = 0

        # Iterate through file line by line.
        for line in fileHandle:

            # First line contains id and name
            if i == 0:
                # must first remove brackets and split on whitespace
                line = line.strip().strip('[').strip(']').split()

                baseElement.id = line[0]
                baseElement.name = line[1]
                    #print(format(i) + ":" + baseElement.id + " " + baseElement.name)    # DEBUG LINE
                # If threshold given, use it for this matrix. Else, use default.
                if (len(line) > 2):
                    baseElement.threshold = float(line[2])
                else:
                    baseElement.threshold = default_th

            # Order of position weight matrices is A,C,G,T. Remove brackets etc.
            elif i < 5:
                # must first remove brackets and split on whitespace
                #     technically should modify index by found ACGT
                    #print(i)    # DEBUG LINE
                line = line.strip().strip('A').strip('C').strip('G').strip('T').strip().strip('[').strip(']').split()
                baseElement.matrix[i - 1] = line[0:]

            # add motif to list and continue (there are 2 newlines between motifs)
            else:
                # check validity and set values passed
                baseElement.checkValid()
                baseElement.base_probability = base_pr
                baseElement.pseudoCount = pc

                if (baseElement.validFlag):
                    # calculate base occurence frequency by position
                    baseElement.calculate_probabilities()

                # QQQ if invalid should it not include?
                # append the current TF matrix to the motif set
                motifSet.addMotif(baseElement)

                # reset the tracker variables
                i = -1
                baseElement.clear()
            i += 1

    return (motifSet)