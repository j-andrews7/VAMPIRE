#!/usr/bin/env python3
"""
sequence.py contains sequence object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality

XXX|YYY: note incomplete
"""


class sequenceArray:
    """this class is the array of sequence Elements"""

    def __init__(self):
        self.seq = []    # array of sequence objects

    def addSeq(self, sequenceElementToAdd):
        """add passed element to the set"""
        self.seq.append(sequenceElementToAdd)
        return

    def deleteSeq(self, seqIndex):
        """
        remove passed element index from the set
        need to handle parentIndex and referenceIndex values
        use linked

        XXX: incomlete/invalid! handling links needs fixed
        """

        if (seqIndex < 0 | seqIndex > len(self.seq)):
            return

        print('WARNING: deleteSeq is incomplete/incorrect. Needs to be fixed')
        return

        handleLinks = False
        if len(self.seq) > 1:
            if len(self.seq.linked) > 0:
                handleLinks = True
        else:
            self.__init__()    # reset to defaults
            return

        # get to here then know started with multi-element set

        # first remove the links
        if handleLinks:
            # for each self.seq.linked
            # must remove link to target to be deleted
            # must shift all index > seqIndex targets down by 1
            # for all parentIndex, referenceIndex, and linked values
            # in all array elements of seq
            print('placeholder function print for index cleanup')

        del self.seq[seqIndex]

        return

    def length(self):
        """return length of the array"""
        return (len(self.seq))

    def setParent(self, childLinking, parentLinked):
        """
        For passed integer indexes link child to parent
        Note: Both parent and child can be snips
        """

        # first check that both are valid
        if (childLinking < 0 | childLinking > len(self.seq)):
            return
        if (parentLinked < 0 | parentLinked > len(self.seq)):
            return

        # -- link them
        # let the parent know the child exists; assumes new ie could be dup
        self.seq[parentLinked].linked.append[childLinking]
        # link child to parent
        self.seq[childLinking].parentIndex = parentLinked

        return

    def setReference(self, childLinking, referenceLinked):
        """For passed integer indexes link 'child' to reference"""

        # first check that both are valid
        if (childLinking < 0 | childLinking > len(self.seq)):
            return
        if (referenceLinked < 0 | referenceLinked > len(self.seq)):
            return

        if (not(self.seq[referenceLinked].reference)):
            print(("Make " + format(referenceLinked) + "reference."))

        # -- link them
        # let the reference know the child exists; assumes new ie could be dup
        self.seq[referenceLinked].linked.append[childLinking]
        # link child to reference
        self.seq[childLinking].referenceIndex = referenceLinked
        # set reference index object as a reference
        self.seq[referenceLinked].reference = True

        return


class sequence:
    """
    class used to define variant and reference sequence objects
    sequence (str): DNA sequence (ACGTN only)

    Some variables(parentIndex, referenceIndex, linked) are use with the
    sequenceArray class

    XXX: change to not use reference linked list
        instead have ref and variant core string + 2 wing strings of max length
        then process each wing and cores separately, combine to compute
        score, this way only score wings once for reference and variant score
        if variant shorter than reference must insert characters from reference
        genome to make the same length
        for each motif just shorten the wings to process

    XXX: also mark the samples for the variant by reading rest of the line
    XXX: build across multiple variant indexes in wings with matching samples

    """

    def __init__(self):
        self.name = ""       # often the chromosome name
        self.position = -1   # position within genome (index)
        self.seq = ""        # sequence as a string
        self.seqInt = []     # sequence as a set of numbers
        self.snip = True     # false if not a sequence snip
        self.positionCenter = -1       # for snips defines snip center point
        self.revComplement = ""        # sequence reverse complement as string
        self.revComplementInt = []     # sequence reverse complement as numbers
        self.parentIndex = -1     # index in array to full sequence
                                  # only applies if snip = True
        self.reference = False    # True if reference genome
        self.referenceIndex = -1  # index in array to reference sequence
                                  # if snip = True, index is to reference snip
                                  # if snip = False, index is to reference
        self.linked = []     # other object indices linked to this object

    def clear(self):
        self = self.__init__()    # ZZZ: note: self = probably not necessary

    def convert2Int(sequence):    # XXX: incomplete
        """
        convert the sequence from character to integer values
        QQQ: what is the motif file order for ACGT? need to match!
        A = 1
        C = 2
        G = 3
        T = 4

        """
        sequenceInt = [] * len(sequence)
        return sequenceInt

    def convert2IntSequence(self):
        """ call convert2Int against the sequence """
        self.seqInt = self.convert2Int(self.seq)
        return

    def convert2IntReverse(self):
        """ call convert2Int against the reverse complement sequence """
        self.revComplementInt = self.convert2Int(self.revComplement)
        return

    def reverse_complement(self):     # QQQ: any reason to not be sequence class method?
        """
        Return reverse complement of a sequence in variable self.revComplement.

        Args:
            sequence (str): DNA sequence (ACGTN only)
            complement: A <--> T, C <--> G

        """
        self.revComplement = ""

        for idx in range(len(self.seq) - 1, -1, -1):
            base = 'N'

            if self.seq[idx].upper() == 'A':
                base = 'T'
            elif self.seq[idx].upper() == 'C':
                base = 'G'
            elif self.seq[idx].upper() == 'G':
                base = 'C'
            elif self.seq[idx].upper() == 'T':
                base = 'A'

            self.revComplement += base

        return

    # def reverse_complementInt     functionality by rev str --> convert2int
    #    i.e. convert2IntReverse()

    def subString(sequenceStr, wing_length):
        """ given a sequenceString return wing_length around positionCenter """
        subString = ""

        return (subString)


# ______ START METHODS RELATED TO THE CLASS BUT NOT TIED TO IT _____

# QQQ: rethink --> need ref and variant sequence for surrounding
def get_surrounding_seq(chromo, var_pos, ref_l, wing_l, fas):    # QQQ: any reason to not be sequence class method?
    """ Return sequence containing variant base + specified number
        of bases on each side from reference sequence file.

    Args:
        chromo = Chromosome of variant, e.g. "chr19" or "19" (must be a string)
        var_pos = Integer position of variant start within chromosome.
        ref_l = length of reference sequence. Will be 1 for SNP/insertions but
            greater than 1 for deletions (e.g. deletion of ACTG to G -> ref_l is 4)
        wing_l = Integer number of bases on each side of variant to return (wing
            length) a s full sequence.
        fas = indexed fasta file object

    Returns:
        ref_seq = Sequence (string) containing the variant base + specified
        number of bases on each side from reference sequence file.
    """

    """#debug
        print("\tSearching indexed fasta chromosome "+str(chromo)+" at pos "+
            str(var_pos - wing_l - 1)+":"+str(var_pos + wing_l + ref_l - 1))"""

    # fas['chr1'][0:1] returns the first base (just one)
    # is either  0 indexed and max (last) base is not returned
    #   or      1 indexed and min (first) base isn't returned
    # A 'sequence' object is returned by fas, but converted by a string for return
    ref_seq = fas[chromo][var_pos - wing_l - 1: var_pos + wing_l + ref_l - 1]

    # debug print("\tSequence: "+str(ref_seq))

    return str(ref_seq)


def readLineToSampleDictionaries(headerString):
    """
    Parse header of VCF file to return dictionary of sample names found in file

    Given a string break on tabs, convert 9th-end elements to two dictionaries
        1) samplesByName:  where word at element is key and value is index
        2) samplesByIndex: where word at element is value and key is index

    Args:
        header_line (str): Header line from VCF file.

    Returns:
        2 dictionaries (see above description)
        (samplesByName, samplesByIndex)

    reference: based on get_vcf_samples
    """

    samplesByName = {}
    samplesByIndex = {}

    line_list = headerString.strip().split("\t")

    samples = line_list[9:]

    for index in range(len(samples)):
        sampleName = samples[index].split(".")[-1]
        if (index > 0) and (sampleName in samplesByName):
            # WARNING: previously defined key Name! QQQ: occurs? If yes --> bad
            print((sampleName + " at " + format(index) + " repeat of index " +
                    samplesByName[sampleName]))
            # make sampleName unique by adding index
            sampleName += format(index)

        # add to the dictionary
        samplesByName[sampleName] = index
        samplesByIndex[index] = sampleName

    return (samplesByName, samplesByIndex)

