#!/usr/bin/env python3

import sequence
import motif
import time
from pyfaidx import Fasta

def timeString():
    """ Return time as a string YearMonthDay.Hour.Minute.Second
    """
    return time.strftime('%Y%m%d.%H:%M:%S')


print(("Importing the motif_set " + timeString()))
file_motif = '../../data/HOCOMOCOv10.JASPAR_FORMAT.TF_IDS.txt'
pc = 0.1
th = 0.225
bp = [0.27, 0.23, 0.28, 0.22]
motif_set = motif.get_motifs(file_motif, pc, th, bp)

print(("Indexing the Reference Genome " + timeString()))
file_reference_genome = '../../data/genome_reference/reference_genome_hg19.fa'
fa_ind = Fasta(file_reference_genome)

print(("Importing the variants " + timeString()))
wing_l = 50
file_input = '../../data/FLDL_CCCB_RARE_VARIANTS.MERGED.RNA_DP10.RNA_NODUPS.CHIP_MULTIMARK.SORTED.vcf'

# read variant file variant set without printing header information
variant_set = sequence.SequenceArray()
with open(file_input) as vcf_handle:

    line = vcf_handle.readline()

    # Skip info lines
    while line.startswith("#"):
        # Print new info lines at the top of the ##INFO section
        line = vcf_handle.readline()

    # Process each variant
    eof_counter = 0
    while True:
        # Reads in the next line of the vcf
        # adds the next variant's information to the element array object

        line = vcf_handle.readline()

        if line is None:    # stop infinite loop
            print(('debug: stop file read with count' + format(eof_counter)))
            break

        # skip empty and information lines
        if line.startswith("##"):
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

print(("Finished Importing the variants " + timeString()))
print(("set size: " + format(len(variant_set.seq))))
print(("set size: " + format(variant_set.length())))
print(("name 0: " + format(variant_set.seq[0].name) + ":" + format(variant_set.seq[0].position)))
print(("name 1: " + format(variant_set.seq[1].name) + ":" + format(variant_set.seq[1].position)))

index = 0
var_element = variant_set.seq[index]
var_element.get_surround_seq(wing_l, fa_ind, False)

t1 = time.time()
var_element.assign_rev_complement()
ref_seq = var_element.return_full_ref_seq(wing_l)
var_seq = var_element.return_full_var_seq(wing_l)
plusmatch = motif_set.motif_match(bp, ref_seq, var_seq, wing_l)
t2 = time.time()
print(("Processing time for original: " + format(t2 - t1)))

return -1    # XXX: incomplete
t1 = time.time()
var_element.assign_int_versions()
ref_seq = var_element.return_full_ref_seq(wing_l)    # need integer version
var_seq = var_element.return_full_var_seq(wing_l)
plusmatch = motif_set.motif_match(bp, ref_seq, var_seq, wing_l)    # need integer version
t2 = time.time()
print(("Processing time for original: " + format(t2 - t1)))





