#!/usr/bin/env python3

import pytest

# this is not a correct test file but worked for interactive testing
import sequence
import sequence
x = sequence.sequence_element()
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

