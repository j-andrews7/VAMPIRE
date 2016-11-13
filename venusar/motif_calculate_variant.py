#!/usr/bin/env python3

import sequence
import motif
import time
from pyfaidx import Fasta
import pdb    # necessary for debugger; use pdb.set_trace()
import numpy as np

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
file_input = '../../data/FLDL_CCCB_RARE_VARIANTS.MERGED.RNA_DP10.RNA_NODUPS.CHIP_MULTIMARK.SORTED-head160.vcf'

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


t3 = time.time()
var_element.assign_int_versions()
ref_seq = var_element.return_full_ref_seq_int(wing_l)    # need integer version
var_seq = var_element.return_full_var_seq_int(wing_l)
plusmatch_int = motif_set.motif_match_int(bp, ref_seq, var_seq, wing_l)    # need integer version
# XXX will need reverse complement versions of integers later
t4 = time.time()
print(("Processing time for int replacement including conversion time: " + format(t4 - t3)))


#import numpy as np
#y = np.array(motif_set.motifs[0].matrix, float)
t5 = time.time()
var_element.assign_int_versions()
ref_seq = var_element.return_full_ref_seq_int(wing_l)    # need integer version
var_seq = var_element.return_full_var_seq_int(wing_l)
plusmatch_np = motif_set.motif_match_np(
    np.array(bp), np.array(ref_seq), np.array(var_seq), wing_l)    # need integer version
# XXX will need reverse complement versions of integers later
t6 = time.time()
print(("Processing time for int-np replacement including conversion time: " + format(t6 - t5)))


print(("comparison ref score [" +
    format(plusmatch[1].ref_score) +
    "] to int [" + format(plusmatch_int[1].ref_score) +
    "] to np [" + format(plusmatch_np[1].ref_score) + "]"
    ))

print(("comparison var score [" +
    format(plusmatch[1].var_score) +
    "] to [" + format(plusmatch_int[1].var_score) + "]"
    "] to np [" + format(plusmatch_np[1].var_score) + "]"
    ))


pdb.set_trace()

# testing
x = np.array([0, 1, 2, 3, 4, 5, 6], int)
sequence.crop_from_left(x, 3)
#array([3, 4, 5, 6])
x2 = sequence.crop_from_left(x, 3)
x2.dtype
#dtype('int64')
y = np.array([0, 1, 2, 3, 4, 5, 6], float)
x * y
#array([  0.,   1.,   4.,   9.,  16.,  25.,  36.])
z = x * y
z.dtype
#dtype('float')
w2 = np.array([0., 1., 2., 3., 4., 5., 6.1], float)
w2.dtype
#dtype('float64')
w2 = np.array([0, 1, 2, 3, 4, 5, 6.1])
w2.dtype
#dtype('float64')
# ref: https://docs.scipy.org/doc/numpy/reference/routines.array-manipulation.html
y = np.array(motif_set.motifs[0].matrix)
y2 = np.vstack((np.repeat(0, y.shape[1]), y))
y2.shape
y.shape
#(4, 9)
y2.shape
#(5, 9)
x = np.array([0, 1, 2, 3, 4, 1])
np.arange(11)
#array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10])
np.arange(1, 11)
#array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10])
np.arange(1, 11.1, dtype=int)
#array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11])
collookup = np.arange(len(x))    # the position
# the x values themselves are the row lookup
y2 = np.asarray(y2 * 100, dtype=int)    # easier to read
y2[x, collookup]
g = motif.MotifElement()
#g.score_base( y2[ x, collookup ], .25 )    # this fails
bp = [0.27, 0.23, 0.28, 0.22]
bp5 = np.insert(np.array(bp), 0, 0)    # with 0 element probability
xbp = bp5[x]                # base probability for ACGTN at position
xvp = y2[x, collookup]    # motif probability for ACGTN at position
#np.apply_along_axis( log2, 0, xbp )    # errors
np.apply_along_axis(motif.log2infforcednp, 0, xvp / xbp)    # could use apply_across_axis with score_base
#array(nan)
# http://stackoverflow.com/questions/21687581/typeerror-only-length-1-arrays-can-be-converted-to-python-scalars-while-trying#21697288


#### why does numpy take so much longer? test the scoring mechanism loop
#    answer: scoring along is twice a slow just for the calculation using numpy!
x = [0, 1, 2, 3, 4, 1]
xnp = np.array(x)
y = [ [.25, .26, .2, .9, .4, .9],
      [.45, .16, .3, .7, .2, .1],
      [.65, .06, .4, .4, .3, .2],
      [.85, .86, .5, .2, .4, .3] ]
ynp = np.array(y)
ynp = np.vstack((np.repeat(0, ynp.shape[1]), ynp))    # with 0 prob row
bp = [0.27, 0.23, 0.28, 0.22]
bpnp = np.insert(np.array(bp), 0, 0)    # with 0 element probability


t7a = time.time()
score = 0
for pos in range(len(x)):
    # Match base: compute score for each overlap position
    ind = x[pos] - 1
    if ind >= 0:
        score += motif.log2infforced(y[ind][pos] / bp[ind])
score
t7b = time.time()
print(("Processing time for int replacement: " + format(t7b - t7a)))


t8a = time.time()
collookup = np.arange(len(x))    # the position
    # collookup uses motif length not len(sequence_np)) should be equal
xbp = bpnp[xnp]            # base probability for ACGTN at position
xvp = ynp[xnp, collookup]  # motif probability for ACGTN at position
scorenp = np.nansum((np.log2(xvp / xbp)))    # this replaces score_base
scorenp
t8b = time.time()
print(("Processing time for int-np replacement: " + format(t8b - t8a)))