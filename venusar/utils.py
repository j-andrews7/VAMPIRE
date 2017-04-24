#!/usr/bin/env python3
"""
utils.py contains object classes and other related methods.
Does NOT contain a main function.

Intended to be imported with motifs.py and other for related functionality.

XXX|YYY: note incomplete.
"""


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
        Return whether self overlaps Position pos_b's wings.

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
