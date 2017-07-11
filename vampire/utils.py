#!/usr/bin/env python3
"""
utils.py contains generic object classes and other related functions.

Intended to be imported to reduce code redundancy and streamline code maintenance.
"""

import time
import sys


class Position(object):
    """
    Use to represent and handle genomic ranges more easily.

    Args:
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
    """

    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start = start_pos
        self.end = end_pos
        self.wing_start = None
        self.wing_end = None

    def overlaps(self, pos_b):
        """
        Return whether self overlaps Position pos_b's start and end.

        Args:
            pos_b (Position): Another Position.

        Returns:
            bool: True if self overlaps with Position pos_b. False if not.
        """
        if pos_b is None:
            return False

        if self.chrom != pos_b.chrom:
            return False

        start1, start2, end1, end2 = (self.start, pos_b.start, self.end, pos_b.end)

        return end1 >= start2 and end2 >= start1

    def overlaps_wings(self, pos_b):
        """
        Return whether self overlaps Position pos_b's start and end.

        Args:
            pos_b (Position): Another Position.

        Returns:
            bool: True if self overlaps with Position pos_b. False if not.
        """
        if pos_b is None:
            return False

        if self.chrom != pos_b.chrom:
            return False

        start1, start2, end1, end2 = (self.start, pos_b.wing_start, self.end, pos_b.wing_end)

        return end1 >= start2 and end2 >= start1

    def set_wings(self, wing_length):
        """
        Set wing positions for a Position object.

        Args:
            position (int): Position to add the wings to.
            wing_length (int): The length of the wings to add to each side of the position.

        Returns:
            wing_positions (tuple): A tuple containing the position of each wing.
        """
        wing_length = int(wing_length)

        self.wing_start = self.start - wing_length
        self.wing_end = self.end + wing_length

        return

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)


def timeString():
    """
    Return time as a string YearMonthDay.Hour.Minute.Second
    """
    return time.strftime('%Y%m%d.%H:%M:%S')


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size
