#!/usr/bin/env python3
# ref: http://docs.python-guide.org/en/latest/writing/tests/

# this is not a correct test file but worked for interactive testing
import pytest
import sequence

# -- get list of functions to test
#    grep -n '^class' sequence.py > sequenceFunctionSet
#    grep -n 'def ' sequence.py >> sequenceFunctionSet
#    cat sequenceFunctionSet | sort -n > f1; mv f1 sequenceFunctionSet
#  12:class SequenceArray:
#  19:    def __init__(self):
#  24:    def add_seq_defined(self, chromosome, position, reference_seq, variant_seq):
#  44:    def add_seq(self, sequence_element_to_add):
#  49:    def add_seq_multivariant(self, sequence_element_to_add):
#  55:    def clear(self):
#  59:    def delete_seq(self, seq_index):
#  84:    def length(self):
#  88:    def multivariant_list_build(self, limit_mv_dist):
#  114:class SequenceElement:
#  132:    def __init__(self):
#  144:    def assign(self, chromosome, position, reference_seq, variant_seq):
#  161:    def assign_rev_complement(self):
#  170:    def assign_samples(self, file_samples_array):
#  192:    def clear(self):
#  195:    def get_surround_seq(self, wing_length, fasta_object, force_ref_match):
#  240:    def print_str(self):
#  251:    def return_full_seq(self, sequence, wing_length):
#  267:    def return_full_seq_reverse_complement(self, sequence, wing_length):
#  284:    def return_full_ref_seq(self, wing_length):
#  291:    def return_full_var_seq(self, wing_length):
#  298:    def return_full_ref_seq_reverse_complement(self, wing_length):
#  305:    def return_full_var_seq_reverse_complement(self, wing_length):
#  313:class SequenceStr:
#x  318:    def __init__(self):
#x 327:    def convert2int(sequence):    # QQQ: convert to numpy array?
#x 353:    def convert2int_sequence(self):
#x 358:    def convert2int_rev_complement(self):
#x 363:    def reverse_complement(self):
#x  390:    # def reverse_complement_int     functionality by rev str --> convert2int
#x 397:def crop_from_left(sequence_string, crop_length):
#x 414:def crop_from_right(sequence_string, crop_length):
#  432:def get_surrounding_seq(chromo, var_pos, ref_l, wing_l, fas):
#  465:def read_line2sample_dictionaries(headerString):
#x 507:def sub_from_left(sequence_string, return_length):
#x 526:def sub_from_right(sequence_string, return_length):

@pytest.fixture
def ss():
    s = sequence.SequenceStr()
    s.seq = "ACGT"
    return s


def ssbad():
    s = sequence.SequenceStr()
    s.seq = "ACGTN"
    return s


def test_convert2int(ss, ssbad):
    ssi = ss.convert2int(ss.seq)
    assert ssi == [1, 2, 3, 4]
    ss.convert2int_sequence()
    assert ssi == ss.seq_int
    ssbad.convert2int_sequence()
    assert ssbad.seq_int == [1, 2, 3, 4, 0]


def test_reverse_complement(ss, ssbad):
    ss.reverse_complement()
    assert ss.seq_rev_complement == 'ACGT'
    ssbad.reverse_complement()
    assert ssbad.seq_rev_complement == 'NACGT'


def test_convert2int_rev_complement(ss, ssbad):
    ss.convert2int_rev_complement()
    assert ss.seq_rev_complement_int == [1, 2, 3, 4]
    ssbad.convert2int_rev_complement()
    assert ssbad.seq_rev_complement_int == [0, 1, 2, 3, 4]

def test_stringfunctions(ss, ssbad):
    w = "Some Sentence Does Not Matter"

    # sub length tests (positive and negative)
    cropn = 5
    assert sequence.crop_from_left(w, cropn) == w[cropn:]
    assert sequence.crop_from_right(w, cropn) == w[:-cropn]
    assert sequence.crop_from_left(w, -cropn) == w
    assert sequence.crop_from_right(w, -cropn) == w

    assert sequence.sub_from_left(w, cropn) == w[-cropn:]
    assert sequence.sub_from_right(w, cropn) == w[:cropn]
    ww = sequence.sub_from_left(w, -cropn)
    assert len(ww) == 0 and ww == ''
    ww = sequence.sub_from_right(w, -cropn)
    assert len(ww) == 0 and ww == ''

    # length tests (positive)
    wlen = len(w)
    ww = sequence.crop_from_left(w, wlen)
    assert ww == w[wlen:] and len(ww) == 0
    ww = sequence.crop_from_right(w, wlen)
    assert ww == w[wlen:] and len(ww) == 0

    ww = sequence.sub_from_left(w, wlen)
    assert ww == w[:wlen] and len(ww) == wlen
    ww = sequence.sub_from_right(w, wlen)
    assert ww == w[:wlen] and len(ww) == wlen

    # greater than length tests (positive)
    cropn = wlen + 1
    assert sequence.crop_from_left(w, cropn) == w[cropn:]
    assert len(sequence.crop_from_right(w, cropn)) == 0
    assert sequence.crop_from_left(w, cropn) == w[cropn:]
    assert len(sequence.crop_from_right(w, cropn)) == 0

    ww = sequence.sub_from_left(w, cropn)
    assert ww == w[:wlen] and len(ww) == wlen
    ww = sequence.sub_from_right(w, cropn)
    assert ww == w[:wlen] and len(ww) == wlen



#  397:def crop_from_left(sequence_string, crop_length):
#  414:def crop_from_right(sequence_string, crop_length):
#  432:def get_surrounding_seq(chromo, var_pos, ref_l, wing_l, fas):
#  465:def read_line2sample_dictionaries(headerString):
#  507:def sub_from_left(sequence_string, return_length):
#  526:def sub_from_right(sequence_string, return_length):



x = sequence.SequenceElement()
x.sub_left_wing("random",3)
#'dom'
x.sub_right_wing("random",3)
#'ran'
x.seq_left_wing.seq = "leftwing"
x.seq_right_wing.seq = "rightwing"
x.return_full_seq( "random", 4 )

# WHY!#$*)(#^#!
#Traceback (most recent call last):
  #File "<stdin>", line 1, in <module>
  #File "/data640g1/data/documents/docs/projects/payton_venusar/VENUSAR_DEV/venusar/sequence.py", line 157, in return_full_seq
    #build_seq = sub_left_wing(self.seq_left_wing.seq, wing_length)
#NameError: name 'sub_left_wing' is not defined
x.clear()
del x
x = sequence.SequenceElement()
x.seq_var.seq = 'var seq'
x.seq_ref.seq = 'ref seq'
x.seq_left_wing.seq = 'left words'
x.seq_right_wing.seq = 'right words'
x.return_full_ref_seq(5)
ref_seq = x.return_full_ref_seq(5)


import sequence
y = sequence.SequenceStr()
y.seq = 'ACGT'
y.convert2int_sequence()
y.seq_int
[1, 2, 3, 4]
y.convert2int(y.seq)
[1, 2, 3, 4]

# -- test SequenceArray capability
yAr = sequence.SequenceArray()
# chromosome, position, reference_seq, variant_seq):
yAr.add_seq_defined('xyz', 32, 'A', 'T')
yAr.add_seq_defined('xyz', 1, 'C', 'G')
yAr.add_seq_defined('xyz', 12, 'C', 'G')
yAr.add_seq_defined('abc', 1000, 'C', 'G')
yAr.add_seq_defined('abc', 100, 'C', 'G')
yAr.add_seq_defined('mn3', 83, 'C', 'G')
yAr.add_seq_defined('abc', 13, 'C', 'G')

print('##  -- pre sort -- ##')
ind = 0
while ind < len(yAr.seq):
    print(((yAr.seq[ind]).name + ":" + format((yAr.seq[ind]).position)))
    ind = ind + 1

yAr.sort()

print('##  -- post sort -- ##')
ind = 0
while ind < len(yAr.seq):
    print(((yAr.seq[ind]).name + ":" + format((yAr.seq[ind]).position)))
    ind = ind + 1


# -- working on multivariant_list_build
x = [0, 5, 3, 4]
y = [1, 3, 5, 6, 8]
overlap = list(set(x).intersection(y))    # returns set; must convert to list
list(set(x).difference(y))                # items in x not in y


file_reference_genome = "../../data/genome_reference/reference_genome_hg19.fa"
from pyfaidx import Fasta
fa_ind = Fasta(file_reference_genome)
# name    position        ref        var
# chr1	14653	.	C	T	48.56	PASS
# chr1	14907	.	A	G	28.47	PASS
# chr1	14930	.	A	G	28.47	PASS
#
start_pos = 14907
next_pos = 14930
chr_name = "chr1"
variant_str = "G"
variant_str2 = "G"
sequence.get_surrounding_seq(chr_name, start_pos, next_pos - start_pos + 1, 0, fa_ind)
# returns: C

ref_seq = sequence.get_surrounding_seq(chr_name, start_pos,
    next_pos - start_pos + len(variant_str2), 0, fa_ind)
#    2. get inbetween ref_seq (excludes variant positions)
between_seq = sequence.get_surrounding_seq(chr_name, start_pos + len(variant_str),
    next_pos - (start_pos + len(variant_str)), 0, fa_ind)
#    3. build the combined variant sequence
var_seq = variant_str + between_seq + variant_str2

ref_seq = get_surrounding_seq(chr_name, start_pos, 1, limit_mv_dist + wing_l, fasta_object)



x = 'AGAAAAAGGCAGGACAGAATTACAN'
import re
search_pattern_notaz = re.compile('[^a-zA-Z]')    # faster to precompile
re.search('[a-zA-Z]', x)    # match alpha
    # <_sre.SRE_Match object; span=(0, 1), match='A'>
re.search('[^a-zA-Z]', x)   # match not alpha
re.search('[0-9]', x)       # match number
re.search('[^0-9]', x)      # match not number
    # <_sre.SRE_Match object; span=(0, 1), match='A'>
re.search('[^a-zA-Z]', x + ",")    # match no alpha
    # <_sre.SRE_Match object; span=(25, 26), match=','>
y = re.search(search_pattern_notaz, x + ",")
y is not None
