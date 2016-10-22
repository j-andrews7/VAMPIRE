#!/usr/bin/env python3
"""
sequence.py contains sequence object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality

XXX|YYY: note incomplete
"""
import re    # regular expression matching

class sequenceArray:
    """
    this class is the array of sequence Elements

    XXX: build across multiple variant indexes in wings with matching samples
    """

    def __init__(self):
        self.seq = []        # array of sequence objects
        self.multivariant = []    # list of multivariant seq indices
        return

    def add_seq_defined(self, chromosome, position, reference_seq, variant_seq):
        """
        Create a sequenceElement object, assign values, then add to the set

        Args:
            chromosome number as a string e.g. "chr1",
            position = position of the variant as an integer
            reference_seq = reference sequence
            variant_seq = variant sequence
        Returns:
            assigns argument inputs to element values, no return updates self

        Note: if using samples likely safer to directly call
            sequenceElement Functions so know correclty grouped samples
        """
        new_sequence_element = sequenceElement()
        new_sequence_element.assign(chromosome, position, reference_seq, variant_seq)
        self.add_seq(new_sequence_element)
        return

    def add_seq(self, sequence_element_to_add):
        """add passed sequence element to the set"""
        self.seq.append(sequence_element_to_add)
        return

    def add_seq_multivariant(self, sequence_element_to_add):
        """add passed multivariant sequence element to the set"""
        self.add_seq(sequence_element_to_add)
        self.multivariant.append(len(self.seq) - 1)
        return

    def clear(self):
        self.__init__()
        return

    def delete_seq(self, seq_index):
        """remove passed element index from the set"""

        if (seq_index < 0 | seq_index > len(self.seq)):
            return

        handleLinks = False
        if len(self.seq) > 1:
            if len(self.multivariant) > 0:
                handleLinks = True
        else:
            self.__init__()    # reset to defaults
            return

        # get to here then know started with multi-element set

        # first remove the multivariant index link
        if handleLinks & seq_index in self.multivariant:
            multi_var_index = self.multivariant.index(seq_index)
            del self.multivariant[multi_var_index]
        # now remove the actual element
        del self.seq[seq_index]

        return

    def length(self):
        """return length of the array"""
        return (len(self.seq))

    def multivariant_list_build(self, limit_mv_dist):
        """
        create a new sequenceElement for joined variant set for each case:
            same chromosome + same sample + within limit_mv_dist --> join seq
        characters between variants drawn from reference genome fasta index

        for the set of sequence elements use position, samples, wing size,
        and user defined overlap extent (limit_mv_dist) integer to join
        multiple chromosome variants for sample that fall within wings

        Search is from low to high position

        returns
            updates the seq array to include multi-variant objects
            add

        """
        print("XXX: function NOT implemented yet. Please check in later.")
        # multivar_str = var_seq1 + reference_space + var_seq2
        # new_element = sequenceElement()
        #    # grab name, position, snip from var_seq1
        # new_element.seq_var = multivar_str
        # sequenceArray.add_seq(new_element)
        return


class sequenceElement:
    """
    class used to define sequence objects built from
    sequenceStr objects for separate variant, reference, and wing strings
    Mark the samples for the variant by reading rest of the line

    How Used:
    Uses ref and variant core string + 2 wing strings of max length
    Does not insert extra empty characters if variant is shorter than reference
    For each motif just shorten the wings to process

    QQQ: Thought process each wing and cores separately, combine to compute
        score, this way only score wings once for reference and variant score

    YYY: note: substrings built so wings have length of motif

    """

    def __init__(self):
        self.name = ""       # often the chromosome name
        self.position = -1   # position within genome (index)
        self.snip = True     # false if not a sequence snip
        self.seq_var = sequenceStr()    # variable object for variant sequence
        self.seq_ref = sequenceStr()    # variable object for reference seq.
        self.seq_left_wing = sequenceStr()   # left wing sequence object
        self.seq_right_wing = sequenceStr()  # right wing sequence object
        self.samples = []     # set of sample (by index) for this sequence
                              # index from placement in file just like header

    def assign(self, chromosome, position, reference_seq, variant_seq):
        """
        Assign values to the element variables

        Args:
            chromosome number as a string e.g. "chr1",
            position = position of the variant as an integer
            reference_seq = reference sequence
            variant_seq = variant sequence
        Returns:
            assigns argument inputs to element values, no return updates self
        """
        self.name = chromosome
        self.position = position
        self.seq_var.seq = reference_seq
        self.seq_ref.seq = variant_seq

    def assign_samples(self, file_samples_array):
        """
        Convert passed sparse array into array of populated items only
        Use in conjunction with read_line2sample_dictionaries()
        Function must assume that array is the same size dictionary sample set

        Args:
            file_samples_array (technically a list)
        Returns:
            Sets self.samples
        """

        #file_samples_array = line_list[9:]    # example of input

        # using regular expression to match on numbers to append sample set
        # QQQ: maybe use lenght not .\. or convert to number or ? maybe faster
        search_pattern = re.compile('[0-9]')    # faster to precompile
        for index in range(len(file_samples_array)):
            if re.search(search_pattern, file_samples_array[index]) is not None:
                self.samples.append(index)

        return

    def clear(self):
        self.__init__()

    def return_full_seq(self, sequence, wing_length):
        """
        return sequence to process with wings attached
        arguments:
            sequence should be the core sequence
            wing_length = integer number of characters to use from the wing seq
        """
        build_seq = self.sub_left_wing(self.seq_left_wing.seq, wing_length)
        build_seq += sequence
        build_seq += self.sub_right_wing(self.seq_right_wing.seq, wing_length)
        return build_seq

    def return_full_seq_reverse_complement(self, sequence, wing_length):
        """
        return sequence to process with wings attached
        arguments:
            sequence should be the reverse complement of the core sequence
            wing_length = integer number of characters to use from the wing seq
        """
        build_seq = self.sub_right_wing(self.seq_right_wing.seq_rev_complement, wing_length)
        build_seq += sequence
        build_seq += self.sub_left_wing(self.seq_left_wing.seq_rev_complement, wing_length)
        return build_seq

    def return_full_ref_seq(self, wing_length):
        """return reference string to process from seq_ref, wings, and lengths"""
        return self.return_full_seq(self.seq_ref.seq, wing_length)

    def return_full_var_seq(self, wing_length):
        """return variant string to process from seq_var, wings, and lengths"""
        return self.return_full_seq(self.seq_var.seq, wing_length)

    def return_full_ref_seq_reverse_complement(self, wing_length):
        """return reference string to process from seq_ref, wings, and lengths"""
        return self.return_full_seq(self.seq_ref.seq_rev_complement, wing_length)

    def return_full_var_seq_reverse_complement(self, wing_length):
        """return variant string to process from seq_var, wings, and lengths"""
        return self.return_full_seq(self.seq_var.seq_rev_complement, wing_length)


    def sub_left_wing(self, sequence_string, return_length):
        """
        return a sequence string of return_length from sequence_string
        drops values at start ie chops off left side
        if return_length > length sequence_string --> does NOT pad
        """
        if (return_length <= 0):
            return ""

        seq_length = len(sequence_string)
        if (return_length > seq_length):
            return sequence_string

        return sequence_string[(seq_length - return_length):]

    def sub_right_wing(self, sequence_string, return_length):
        """
        return a sequence string of return_length from sequence_string
        drops values at end ie chops off right side
        if return_length > length sequence_string --> does NOT pad
        """
        if (return_length <= 0):
            return ""

        seq_length = len(sequence_string)
        if (return_length > seq_length):
            return sequence_string

        return sequence_string[:return_length]

    def set_seq_left_wing(self, wing_length, fasta_index):
        """
        uses variant position to pull left wing from reference genome
        Args:
            wing_length    integer of number of bases to return
            fasta_index    the fasta indexed reference genome
        """
        # XXX: incomplete
        # self.seq_left_wing.seq = get_surrounding_seq ...
        return


class sequenceStr:
    """
    create sequence string data structure used to build sequence elements
    all sequence strings are DNA sequence (ACGTN only)
    """
    def __init__(self):
        # create the object: number in comment corresponds to many
        #    class function arguments: str_index
        self.seq = ""                   # 0. sequence as a string
        self.seq_int = []               # 1. sequence as a set of numbers
        self.seq_rev_complement = ""    # 2. sequence reverse complement as string
        self.seq_rev_complement_int = []  # 3. sequence reverse complement as numbers
        # QQQ: make numpy arrays of Int?

    def convert2int(sequence):    # QQQ: convert to numpy array?
        """
        convert the sequence from character to integer values
            A = 1
            C = 2
            G = 3
            T = 4
            0 used for all other values
        Args: string
        Returns: array of numbers same length as string
        """
        sequence_int = [] * len(sequence)
        for idx in range(len(sequence)):
            if sequence[idx].upper() == 'A':
                sequence_int[idx] = 1
            elif sequence[idx].upper() == 'C':
                sequence_int[idx] = 2
            elif sequence[idx].upper() == 'G':
                sequence_int[idx] = 3
            elif sequence[idx].upper() == 'T':
                sequence_int[idx] = 4
            else:
                sequence_int[idx] = 0

        return sequence_int

    def convert2int_sequence(self):
        """ call convert2int against the sequence """
        self.seq_int = self.convert2int(self.seq)
        return

    def convert2int_rev_complement(self):
        """ call convert2Int against the reverse complement sequence """
        self.seq_rev_complement_int = self.convert2int(self.seq_rev_complement)
        return

    def reverse_complement(self):
        """
        Set reverse complement of sequence as variable self.seq_rev_complement.

        Args:
            self.seq is sequence (str): DNA sequence (ACGTN only)
            complement: A <--> T, C <--> G

        """
        self.seq_rev_complement = ""

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

            self.seq_rev_complement += base

        return

    # def reverse_complement_int     functionality by rev str --> convert2int
    #    i.e. convert2IntReverse()


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


def read_line2sample_dictionaries(headerString):
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
    QQQ: may be better to use list and use var.index(name) to get index
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

