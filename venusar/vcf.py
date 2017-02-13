#!/usr/bin/env python3
"""
vcf.py contains vcf file manipulation routines
mostly focused on reading

does NOT contain a main function
does NOT contain a class because does NOT currently create a VCF object

intended to be imported and used by other code in venusar package
"""

import sequence


def read_vcf_data(vcf_file, samples_as_dictionary=False):
    """
    Read vcf file, return samples from header, a variant SequenceArray and ...

    Args:
        vcf_file (str): valid vcf file location
        samples_as_dictionary (boolean, False): if True return samples as list
            if false return a list of 2 dictionaries: (samplesByName, samplesByIndex)
            True: calls sequence.read_line2sample_list
            False: calls sequence.read_line2sample_dictionaries
    Returns:
        variant_set = SequenceArray object of the variants
        samples (list): type modified by samples_as_dictionary boolean
    """

    variant_set = sequence.SequenceArray()
    vcf_samples = []

    bad_file = False
    # Open and Read VCF file: populates SequenceArray, assumes set fits in memory
    with open(vcf_file) as vcf_handle:
        line = vcf_handle.readline().strip()

        # Skip info lines
        while line.startswith("##"):
            line = vcf_handle.readline()

        # Parse VCF sample header line to get samples present in file.
        if line.startswith("#"):
            line = line.strip()
            if samples_as_dictionary:
                #(samplesByName, samplesByIndex) = sequence.read_line2sample_dictionaries(line)
                vcf_samples = sequence.read_line2sample_dictionaries(line)
            else:
                vcf_samples = sequence.read_line2sample_list(line)
        else:
            print(("ERROR: read_vcf_data(): Early Exit due to Incorrect file format.\n\tFile(" +
                format(vcf_file) + ") does not have # line following ##info lines."))
            bad_file = True    # flag to not process rest of the file; effectively returns

        if not bad_file:
            variant_set = read_vcf_variant_lines(vcf_handle)

    return (variant_set, vcf_samples)


def read_vcf_variant_lines(file_handle, verbose_flag=False):
    """
    Given an open file handle object (file_handle), read the variant lines and
    return a sequence.SequenceArray() object

    Args:
        file_handle    open file object
        verbose_flag (boolean, False): if True then prints addition statements
    Return:
        sequence.SequenceArray() object
    """

    variant_set = sequence.SequenceArray()

    # Process each variant; reads the variant lines in a loop
    eof_counter = 0
    while True:
        # Reads in the next line of the vcf
        # adds the next variant's information to the element array object

        line = file_handle.readline()

        if line is None:    # stop infinite loop
            print(('debug: stop file read with count' + format(eof_counter)))
            break

        # skip empty, information, and comment lines
        if line.startswith("#"):
            continue
        if line == "":
            # handle while True infinite loop try to find End of File
            if (eof_counter > 1):
                break
            else:
                eof_counter += 1

        line = line.strip()

        # skip empty lines (after removing arbitrary whitespace characters)
        if line == "":
            continue

        # -- process valid variant lines
        line_list = line.split("\t")
        # variant creation: add_seq_defined(chromosome, position, reference_seq, variant_seq):
        # variant_set.add_seq_defined(line_list[0], int(line_list[1]), line_list[3], line_list[4])
        new_sequence_element = sequence.SequenceElement()
        new_sequence_element.assign(line_list[0], int(line_list[1]), line_list[3], line_list[4])
        # grab samples for variant
        new_sequence_element.assign_samples(line_list[9:])
        # push full line (memory hog, but allows multivar computation
        #    w/o significant code manipulation to track which line is current
        #    already processed, etc.)
        new_sequence_element.vcf_line = line
        variant_set.add_seq(new_sequence_element)
        if verbose_flag:
            print(("\tadd element " + new_sequence_element.name))

    return (variant_set)


# XXX add code to read MOTIF Scores from the 7th tab separated element; see update_vcf_motifs_info
# QQQ is the motifs.py print_peak code ever hit? check output
# QQQ would it be better to use activity.py variant class to accomplish this objective?
#    or modify activity.py to use the full sequence class objects and vcf functions?
# XXX add code to read activity.py:get_variant_output() information; SAMPS* and LOCI* with zscores