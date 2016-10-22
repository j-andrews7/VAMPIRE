#!/usr/bin/env python3

import pytest

# this is not a correct test file but worked for interactive testing
import motif

y = motif.motifArray()
for i in range(4):
    x = motif.motifElement()
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

x = motif.motifElement()
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
