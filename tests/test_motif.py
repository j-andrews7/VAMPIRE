#!/usr/bin/env python3

# this is not a correct test file but worked for interactive testing
import pytest
import motif

# -- get list of functions to test
#    grep -n '^class' motif.py > motifFunctionSet
#    grep -n 'def ' motif.py >> motifFunctionSet
#    cat motifFunctionSet | sort -n > f1; mv f1 motifFunctionSet
# function set
#  14:class MotifArray:
#  17:    def __init__(self):
#  23:    def add_motif(self, motif_element2add):
#  32:    def delete_motif(self, motif_index):
#  53:    def join(self, motif_array2):
#  80:    def length(self):
#  84:    def motif_match(self, baseline_p, ref_seq, var_seq, wing_l):
#  136:    def motif_scores(self, baseline_p, sequence_str, wing_l, normalize):
#  207:    def process_local_env(self, baseline_p, matches, seq_element, co_binders_dict, v_seq, r_seq):
#  310:    def search_set_max_position(self):
#  327:class MotifElement:
#  334:    def __init__(self):
#  350:    def calculate_probabilities(self):
#  394:    def check_valid(self):
#  411:    def clear(self):
#  415:    def ht_matches_in(self, baseline_p, sequence_str):
#  436:    def print_str(self):
#  449:    def score_base(self, probability, baseline_p):
#  462:    def score_motif(self, baseline_p, sequence_str):
#  506:class MotifMatch:
#  513:    def __init__(self, name, var_score, ref_score, motif_index):
#  541:def get_motifs(motif_filename, pc, default_th, base_pr):
#  629:def get_baseline_probs(baseline_f):


@pytest.fixture
def bp():
    # test file read separately
    # just return list
    return [.25, .25, .25, .25]


def matrix(size):
    m = [[] for x in range(4)]
    s = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    if size > len(s):
        size = len(s)
    m[0] = s[:size]
    m[1] = s[:size]
    m[2] = s[:size]
    m[3] = s[:size]


def motif_element():
    me = motif.MotifElement()
    me.name = "chr37"
    me.matrix = matrix(5)
    return


def test_checkvalid(motif_element):
    assert motif_element.check_valid() is True
    assert motif_element.positions == 5


y = motif.MotifArray()
for i in range(4):
    x = motif.MotifElement()
    x.positions = i
    y.add_motif(x)

y.motifs[3].positions
y.max_positions
y.max_positions_index

x.positions = 13
y.add_motif(x)
y.max_positions
y.max_positions_index

x.positions = 2
y.add_motif(x)
y.max_positions
y.max_positions_index

y.length()

y.delete_motif(4)
y.max_positions
y.max_positions_index
y.length()

x = motif.MotifElement()
x.name = 'test'
x.positions = 32

print((x.name + ":" + format(x.positions)))
x.clear()
print((x.name + ":" + format(x.positions) + " post clear"))

#x.matrix = [[] for x in range(4)]
x.valid_flag
x.check_valid()
x.valid_flag
x.matrix[1] = [3, 4, 5]
x.matrix[2] = [3, 4, 5, 6]
x.check_valid()
x.valid_flag


# -- must create files in path using touch, and text editor for test to work
import motif
x = motif.get_baseline_probs( 'emptyfile' )
#**ERROR** Empty file found.
#File should contain only [ PrA PrC PrG PrT ]
#Where PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)
#Continuing with default probabilities:[0.25, 0.25, 0.25, 0.25]

x = motif.get_baseline_probs( 'notemptyfile' )
#**ERROR** Baseline probability file incorrectly formatted.
#File should contain only [ PrA PrC PrG PrT ]
#Where PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)
#Continuing with default probabilities:[0.25, 0.25, 0.25, 0.25]
x = motif.get_baseline_probs( 'notemptyfile_correct')
