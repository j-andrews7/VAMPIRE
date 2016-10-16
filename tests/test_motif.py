#!/usr/bin/env python3

import pytest

# this is not a correct test file but worked for interactive testing
import motif

y = motif.motifArray()
for i in range(4):
    x = motif.motifElement()
    x.positions = i
    y.addMotif(x)

y.motifs[3].positions
y.max_positions
y.max_positionsIndex

x.positions = 13
y.addMotif(x)
y.max_positions
y.max_positionsIndex

x.positions = 2
y.addMotif(x)
y.max_positions
y.max_positionsIndex

y.length()

y.deleteMotif(4)
y.max_positions
y.max_positionsIndex
y.length()

x = motif.motifElement()
x.name = 'test'
x.positions = 32

print((x.name + ":" + format(x.positions)))
x.clear()
print((x.name + ":" + format(x.positions) + " post clear"))

#x.matrix = [[] for x in range(4)]
x.validFlag
x.checkValid()
x.validFlag
x.matrix[1] = [3, 4, 5]
x.matrix[2] = [3, 4, 5, 6]
x.checkValid()
x.validFlag

