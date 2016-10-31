#!/usr/bin/env python3
"""
sequence.py contains sequence object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality

XXX|YYY: note incomplete
"""
import re    # regular expression matching

class SequenceArray:
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
        Create a SequenceElement object, assign values, then add to the set

        Args:
            chromosome number as a string e.g. "chr1",
            position = position of the variant as an integer
            reference_seq = reference sequence
            variant_seq = variant sequence
        Returns:
            assigns argument inputs to element values, no return updates self

        Note: if using samples likely safer to directly call
            SequenceElement Functions so know correclty grouped samples
        """
        new_sequence_element = SequenceElement()
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
        create a new SequenceElement for joined variant set for each case:
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
        # new_element = SequenceElement()
        #    # grab name, position, snip from var_seq1
        # new_element.seq_var = multivar_str
        # SequenceArray.add_seq(new_element)
        return


class SequenceElement:
    """
    class used to define sequence objects built from
    SequenceStr objects for separate variant, reference, and wing strings
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
        self.seq_var = SequenceStr()    # variable object for variant sequence
        self.seq_ref = SequenceStr()    # variable object for reference seq.
        self.seq_left_wing = SequenceStr()   # left wing sequence object
        self.seq_right_wing = SequenceStr()  # right wing sequence object
        self.samples = []     # set of sample (by index) for this sequence
                              # index from placement in file just like header
        self.vcf_line = ""    # vcf input line for variant

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
        self.seq_var.seq = variant_seq
        self.seq_ref.seq = reference_seq

    def assign_int_versions(self):
        """
        call convert2int to compute the integer version of the 4 elements
        (var, ref, left, right) and their reverse complements
        assigns values.

        WARNING: Assumes self.assign_rev_complement() was already called!
        """
        # convert the sequences to integers
        self.seq_var.convert2int()
        self.seq_ref.convert2int()
        self.seq_left_wing.convert2int()
        self.seq_right_wing.convert2int()
        # convert the reverse complements to integers
        self.seq_var.convert2int_rev_complement()
        self.seq_ref.convert2int_rev_complement()
        self.seq_left_wing.convert2int_rev_complement()
        self.seq_right_wing.convert2int_rev_complement()

    def assign_rev_complement(self):
        """
        compute the reverse complement of the 4 elements (var, ref, left, right)
        """
        self.seq_var.reverse_complement()
        self.seq_ref.reverse_complement()
        self.seq_left_wing.reverse_complement()
        self.seq_right_wing.reverse_complement()

    def assign_samples(self, file_samples_array):
        """
        Convert passed sparse array into array of populated items only
        Use in conjunction with read_line2sample_dictionaries()
        Function must assume that array is the same size dictionary sample set

        Args:
            file_samples_array (technically a list)
            file_samples_array = line_list[9:]    # example input from vcf file
        Returns:
            Sets self.samples
        """

        # using regular expression to match on numbers to append sample set
        # QQQ: maybe use length not .\. or convert to number or ? maybe faster
        # YYY: Maybe, so long as the output isn't changed. Could also do:
        #    `any(char.isdigit() for char in file_samples_array[index])` if you
        #     don't want to precompile a regex for every call. No idea on speed.
        # CCC-WK: input file read time with this code was well under a minute,
        #     non-issue, leave for future possible speed increases. XXX double check it was used
        search_pattern = re.compile('[0-9]')    # faster to precompile
        for index in range(len(file_samples_array)):
            if re.search(search_pattern, file_samples_array[index]) is not None:
                self.samples.append(index)

        return

    def clear(self):
        self.__init__()

    def get_surround_seq(self, wing_length, fasta_object, force_ref_match):
        """
        Use general get_surrounding_seq function to query reference genome file
        for the surrounding sequence. Then pull out substrings to define the
        left and right wings. Read then splits because file read is expensive.
        Need wing_length and FASTA index. Other inputs defined.

        Args:
            wing_length = Integer number of bases on each side of variant to
                return full sequence overlapping the variant.
            fasta_object = indexed fasta file object
            force_ref_match = (boolean) if true then force variant reference
                bases to match FASTA reference
        WARNING: if force_ref_match = True and this function was called after
            manipulation of non-seq parts of seq_ref SequenceStr object then
            those parts (int, rev_complement, rev_complement_int) are invalid.
            This function does not handle updating none string sequence parts.
        Returns:
            sets self.seq_left_wing and self.seq_right_wing
        """
        #ref_seq = get_surrounding_seq(chr, pos, len(ref_bases), wing_l, fa_ind)
        surr_seq = get_surrounding_seq(self.name, self.position,
                        len(self.seq_ref.seq), wing_length, fasta_object)

        # check that reference sequence matches vcf's sequence for position
        returned_ref_b = surr_seq[wing_length:-wing_length]
        #print(("sequence: " + surr_seq + "\n\t" + returned_ref_b + "\n" +
        #    "\tleft:  " + sub_from_start(surr_seq, wing_length) + "\n" +
        #    "\tright: " + sub_from_end(surr_seq, wing_length) ))
        if self.seq_ref.seq.upper() != returned_ref_b.upper():
            print(("**ERROR**\nVCF reference sequence for " +
                self.name + " pos " + str(self.position) + ":\n\t\"" +
                self.seq_ref.seq +
                "\" does not match reference file sequence:\n\t\"" +
                returned_ref_b + "\"."))
            if force_ref_match:
                self.seq_ref.seq = returned_ref_b
                print(("Forced reference sequence match for variant " +
                    self.name + " pos " + str(self.position)))

        # assign wing sequence strings
        self.seq_left_wing.seq = sub_from_start(surr_seq, wing_length)
        self.seq_right_wing.seq = sub_from_end(surr_seq, wing_length)
        return

    def print_str(self):
        """
        create a string of the sequence strings to print
        QQQ: add an integer string version later by boolean flag?
        """
        print_string = "\tvariant core   : " + self.seq_var.seq + "\n"
        print_string += "\treference core: " + self.seq_ref.seq + "\n"
        print_string += "\tleft_wing : " + self.seq_left_wing.seq + "\n"
        print_string += "\tright_wing: " + self.seq_right_wing.seq + "\n"
        return print_string

    def return_full_seq(self, sequence_str, wing_length):
        """
        return sequence to process with wings attached. Uses sub_from_end and
        sub_from_start methods to reduce wing sizes to requested length.

        Arguments:
            sequence should be the core sequence
            wing_length = integer number of characters to use from the wing seq
        Returns:
            string of left wing + sequence string + right wings
        """
        build_seq = sub_from_end(self.seq_left_wing.seq, wing_length)
        build_seq += sequence_str
        build_seq += sub_from_start(self.seq_right_wing.seq, wing_length)
        return build_seq

    def return_full_seq_reverse_complement(self, sequence_str, wing_length):
        """
        return sequence to process with wings attached. Uses sub_from_end and
        sub_from_start methods to reduce wing sizes to requested length.

        Arguments:
            sequence should be the reverse complement of the core sequence
            wing_length = integer number of characters to use from the wing seq
        Returns:
            For all strings the reverse complement returns
            string of right wing + sequence string + left wings
        """
        build_seq = sub_from_start(self.seq_right_wing.seq_rev_complement, wing_length)
        build_seq += sequence_str
        build_seq += sub_from_end(self.seq_left_wing.seq_rev_complement, wing_length)
        return build_seq

    def return_full_ref_seq(self, wing_length):
        """
        return reference string to process from seq_ref and wings cropped to
        the specified length
        """
        return self.return_full_seq(self.seq_ref.seq, wing_length)

    def return_full_var_seq(self, wing_length):
        """
        return variant string to process from seq_var and wings cropped to
        the specified length
        """
        return self.return_full_seq(self.seq_var.seq, wing_length)

    def return_full_ref_seq_reverse_complement(self, wing_length):
        """
        return the reverse complement of the reference string to process
        built from from seq_ref and wings cropped to specified length
        """
        return self.return_full_seq(self.seq_ref.seq_rev_complement, wing_length)

    def return_full_var_seq_reverse_complement(self, wing_length):
        """
        return the reverse complement of the variant string to process
        built from from seq_var and wings cropped to specified length
        """
        return self.return_full_seq(self.seq_var.seq_rev_complement, wing_length)


class SequenceStr:
    """
    create sequence string data structure used to build sequence elements
    all sequence strings are DNA sequence (ACGTN only)
    """
    def __init__(self):
        # create the object: number in comment corresponds to many
        #    class function arguments: str_index
        self.seq = ""                   # sequence as a string
        self.seq_int = []               # sequence as a set of numbers
        self.seq_rev_complement = ""    # sequence reverse complement as string
        self.seq_rev_complement_int = []  # sequence reverse complement as numbers
        # QQQ: make numpy arrays of Int?

    def convert2int(self, sequence_str):    # QQQ: convert to numpy array?
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
        sequence_int = [0] * len(sequence_str)
        for idx in range(len(sequence_str)):
            if sequence_str[idx].upper() == 'A':
                sequence_int[idx] = 1
            elif sequence_str[idx].upper() == 'C':
                sequence_int[idx] = 2
            elif sequence_str[idx].upper() == 'G':
                sequence_int[idx] = 3
            elif sequence_str[idx].upper() == 'T':
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


# ______ START METHODS RELATED TO THE CLASS BUT NOT EXPLICITLY TIED TO IT _____


def crop_from_left(sequence_string, crop_length):
    """
    return a sequence string without first crop_length characters
    drops values at start ie chops off left side
    if crop_length > length sequence_string --> returns ""
    if crop_length <= 0 returns string
    """
    if (crop_length <= 0):
        return sequence_string

    seq_length = len(sequence_string)
    if (crop_length > seq_length):
        return ""

    return sequence_string[crop_length:]


def crop_from_right(sequence_string, crop_length):
    """
    return a sequence string without last crop_length characters
    drops values at end ie chops off right side
    if crop_length > length sequence_string --> returns ""
    if crop_length <= 0 returns string
    """
    if (crop_length <= 0):
        return sequence_string

    seq_length = len(sequence_string)
    if (crop_length > seq_length):
        return ""

    return sequence_string[:-crop_length]


# ZZZ: any reason to not be SequenceElement class method?
# YYY: Seems a natural place for it to go to me.
# CCC-WK: leaving outside because does not require sequenceElement parts
def get_surrounding_seq(chromo, var_pos, ref_l, wing_l, fas):
    """ Return sequence containing variant base + specified number
        of bases on each side from reference sequence file.

    Args:
        chromo = Chromosome of variant, e.g. "chr19" or "19" (must be a string)
        var_pos = Integer position of variant start within chromosome.
        ref_l = length of reference sequence. Will be 1 for SNP/insertions but
            greater than 1 for deletions (e.g. deletion of ACTG to G -> ref_l is 4)
        wing_l = Integer number of bases on each side of variant to return (wing
            length) as full sequence.
        fas = indexed fasta file object

    Returns:
        ref_seq = Sequence (string) containing the variant base + specified
        number of bases on each side of variant. Data from reference genome file.
    """

    """#debug
        print("\tSearching indexed fasta chromosome "+str(chromo)+" at pos "+
            str(var_pos - wing_l - 1)+":"+str(var_pos + wing_l + ref_l - 1))"""

    # fas['chr1'][0:1] returns the first base (just one)
    # is either  0 indexed and max (last) base is not returned
    #   or      1 indexed and min (first) base isn't returned
    # fas returns a 'sequence' object, converted by a string for return
    print(("pulling reference sequence for (" + chromo + ", " +
        format(var_pos) + ", " + format(ref_l) + ", " + format(wing_l) + ", " +
        format(fas) + ")\n\t" +
        format(var_pos - wing_l - 1) + ":" + format(var_pos + wing_l + ref_l - 1)
        ))
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
    CCC-JA: This could be an issue if multiple columns are for the same sample.
        This might occur if a user is using multiple methods to call their variants
        and pooling the end results together (which is pretty common and something
        we do as well.) Or they may be calling variants from multiple files for the
        same sample (think replicates or sequencing from different experiments).
        Really, it doesn't matter, so long as we know which samples have the variant
        and which don't. We'll have to create a specific format for the sample names
        in the header so that users can still include extra info if they want (as you
        can see in our input VCF).
    CCC-WK: having non-unique columns is a large problem.
        it will break this dictionary build method as you noted
        it will break lists in function
        QQQ: append .1, .2, .N to each column header to make names unique?
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


def sub_from_end(sequence_string, return_length):
    """
    return a sequence string of return_length from sequence_string
    drops values at start ie chops off left side
    if return_length > length sequence_string --> does NOT pad

    note: function was originally part of SequenceElement as sub_left_wing

    Example:
        sub_from_end('1234567', 3) returns 567
    """
    if (return_length <= 0):
        return ""

    seq_length = len(sequence_string)
    if (return_length > seq_length):
        return sequence_string

    return sequence_string[(seq_length - return_length):]


def sub_from_start(sequence_string, return_length):
    """
    return a sequence string of return_length from sequence_string
    drops values at end ie chops off right side
    if return_length > length sequence_string --> does NOT pad

    note: function was originally part of SequenceElement as sub_right_wing

    Example:
        sub_from_start('1234567', 3) returns 123
    """
    if (return_length <= 0):
        return ""

    seq_length = len(sequence_string)
    if (return_length > seq_length):
        return sequence_string

    return sequence_string[:return_length]

