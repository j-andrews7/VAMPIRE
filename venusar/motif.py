#!/usr/bin/env python3
"""
motif.py contains motif object class and other related methods
does NOT contain a main function

intended to be imported with motifs.py and other for related functionality
"""

import sequence    # need sub_from_*, crop_from_*
                   # (maybe later need int func and split wing/core processing)
                   # then would want base sequence object
from math import log2
import numpy as np

class MotifArray:
    """this class is the array of motif Elements"""

    def __init__(self):
        self.max_positions = 0    # largest positions value for all Elements
        self.max_positions_index = -1    # index of maximum position length
        self.motifs = []
        return

    def add_motif(self, motif_element2add):
        """add passed element to the set"""
        self.motifs.append(motif_element2add)
        if (motif_element2add.positions > self.max_positions):
            self.max_positions = motif_element2add.positions
            self.max_positions_index = len(self.motifs) - 1

        return

    def delete_motif(self, motif_index):
        """remove passed element index from the set"""
        if (motif_index < 0 | motif_index > len(self.motifs)):
            return

        reset_max_length = False
        if len(self.motifs) > 1:
            if motif_index == self.max_positions_index:
                reset_max_length = True
        else:
            self.__init__()    # reset to defaults
            return

        # get to here then know started with multi-element set
        del self.motifs[motif_index]

        if reset_max_length:
            self.search_set_max_position()

        return

    def element_positions_list(self, return_string):
        """
        debug function; return list of MotifArray motifElement position sizes
        Args:
            return_string (boolean): true --> return string; false --> return list
        """
        if return_string:
            return_var = ""
        else:
            return_var = []

        for me in self.motifs:
            if return_string:
                return_var = return_var + "," + format(me.positions)
            else:
                return_var.append(me.positions)

        if return_string:
            return_var.lstrip(",")    # drop leading ,

        return return_var

    def join(self, motif_array2):
        """join another MotifArray object to this MotifArray object"""

        # check validity of the items to be joined
        if type(motif_array2) is not MotifArray:
            return
        if motif_array2.length() == 0:
            return

        # empty current, non-empty new --> replace
        if self.length() == 0:
            self = motif_array2    # allowed in python?
            return

        # comparison elements
        current_max = [self.max_positions, self.max_positions_index]
        joined_max = [motif_array2.max_positions, motif_array2.max_positions_index]
        original_length = self.length()

        # -- compare and join
        self.motifs.append(motif_array2.motifs)
        if (joined_max[0] > current_max[0]):
            self.max_positions = joined_max[0]
            self.max_positions_index = original_length + joined_max[1]

        return

    def length(self):
        """return length of the set of motifs for the array"""
        return (len(self.motifs))

    def motif_match(self, baseline_p, ref_seq, var_seq, wing_l):
        """
        Takes a reference and variant sequence string, then checks for matches
        in the motif set. Outputs match score for both ref seq and var seq
        for any case where either matches above the motif element threshold.
        (former name: match_motifs, wanted grouped in class function list)

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            ref_seq = Reference sequence represented as a string (ACGT only)
            var_seq = Variant sequence represented as a string (ACGT only)
            wing_l = Integer length of sequence of bases flanking the variant
                (generally >= self.max_positions, should match value used to
                create ref_seq and var_seq)

        Returns: A list of MotifMatch objects

        see also: motif_match_int
        """

        # scoring returns: motif match array with the form
        #    (id, name, threshold, max_score, match_seq)
        scored = self.motif_scores(baseline_p, var_seq, wing_l, True)
        r_scored = self.motif_scores(baseline_p, ref_seq, wing_l, True)

        matches = []    # list of MotifMatch objects

        # generate list of motifs that matched var seq to compare ref seq matches
        for idx in range(len(scored)):
            (iden, name, th, max_score, match_seq, motif_index) = scored[idx]
            (rid, rname, rth, rmax_score, rmatch_seq, motif_index) = r_scored[idx]
            if (max_score >= th or rmax_score >= rth):
                if (iden != rid or name != rname or th != rth):
                    print(("***ERROR*** matching motifs to varseq and refseq desynced\n" +
                          iden + " != " + rid +
                          " or " + name + " != " + rname +
                          " or " + format(th) + " != " + format(rth)))
                # tup = (id, name, max_score, match_seq, rmax_score, rmatch_seq)
                match = MotifMatch(name, max_score, rmax_score, motif_index)
                matches.append(match)

        # output variant seq matches vs ref seq matches
        """for m in matches:
            line = "Variant matched motif "+m[1]+" ("+m[0]+") with score "
            line+= str(m[2])[:6]+" at seq: "+m[3]+"\nReference matched motif "+m[1]
            line+= "           with score "+str(m[4])[:6]+" at seq: "+m[5]
            #debug line+= "\n\tRefseq: "+ref_seq+"\n\tVarseq: "+var_seq
            print(line,file=output_f)"""    # WARNING: debug using undefined global

        return matches

    def motif_match_int(self, baseline_p, ref_seq, var_seq, wing_l):
        """
        Takes a reference and variant sequence integer representation,
        then checks for matches in the motif set.
        Outputs match score for both ref seq and var seq
        for any case where either matches above the motif element threshold.
        (former name: match_motifs, wanted grouped in class function list)

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            ref_seq = Reference sequence integer representation (0-4 only)
            var_seq = Variant sequence integer representation (0-4 only)
            wing_l = Integer length of sequence of bases flanking the variant
                (generally >= self.max_positions, should match value used to
                create ref_seq and var_seq)

        Returns: A list of MotifMatch objects

        see also: motif_match
        """

        # scoring returns: motif match array with the form
        #    (id, name, threshold, max_score, match_seq)
        scored = self.motif_scores_int(baseline_p, var_seq, wing_l, True)
        r_scored = self.motif_scores_int(baseline_p, ref_seq, wing_l, True)

        matches = []    # list of MotifMatch objects

        # generate list of motifs that matched var seq to compare ref seq matches
        for idx in range(len(scored)):
            (iden, name, th, max_score, match_seq, motif_index) = scored[idx]
            (rid, rname, rth, rmax_score, rmatch_seq, motif_index) = r_scored[idx]
            if (max_score >= th or rmax_score >= rth):
                if (iden != rid or name != rname or th != rth):
                    print(("***ERROR*** matching motifs to varseq and refseq desynced\n" +
                          iden + " != " + rid +
                          " or " + name + " != " + rname +
                          " or " + format(th) + " != " + format(rth)))
                # tup = (id, name, max_score, match_seq, rmax_score, rmatch_seq)
                match = MotifMatch(name, max_score, rmax_score, motif_index)
                matches.append(match)

        # output variant seq matches vs ref seq matches
        """for m in matches:
            line = "Variant matched motif "+m[1]+" ("+m[0]+") with score "
            line+= str(m[2])[:6]+" at seq: "+m[3]+"\nReference matched motif "+m[1]
            line+= "           with score "+str(m[4])[:6]+" at seq: "+m[5]
            #debug line+= "\n\tRefseq: "+ref_seq+"\n\tVarseq: "+var_seq
            print(line,file=output_f)"""    # WARNING: debug using undefined global

        return matches

    def motif_match_np(self, baseline_p, ref_seq, var_seq, wing_l):
        """
        Takes a reference and variant sequence integer representation,
        then checks for matches in the motif set. ** This version uses numpy **
        Outputs match score for both ref seq and var seq
        for any case where either matches above the motif element threshold.
        This version uses integer representations of the baseline_p and
        sequences passed as numpy arrays.
        (former name: match_motifs, wanted grouped in class function list)

        Args:
            baseline_p = np array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            ref_seq = Reference sequence integer representation (0-4 only)
            var_seq = Variant sequence integer representation (0-4 only)
            wing_l = Integer length of sequence of bases flanking the variant
                (generally >= self.max_positions, should match value used to
                create ref_seq and var_seq)

        Returns: A list of MotifMatch objects
        """

        # scoring returns: motif match array with the form
        #    (id, name, threshold, max_score, match_seq)
        scored = self.motif_scores_np(baseline_p, var_seq, wing_l, True)
        r_scored = self.motif_scores_np(baseline_p, ref_seq, wing_l, True)

        matches = []    # list of MotifMatch objects

        # generate list of motifs that matched var seq to compare ref seq matches
        for idx in range(len(scored)):
            (iden, name, th, max_score, match_seq, motif_index) = scored[idx]
            (rid, rname, rth, rmax_score, rmatch_seq, motif_index) = r_scored[idx]
            if (max_score >= th or rmax_score >= rth):
                if (iden != rid or name != rname or th != rth):
                    print(("***ERROR*** matching motifs to varseq and refseq desynced\n" +
                          iden + " != " + rid +
                          " or " + name + " != " + rname +
                          " or " + format(th) + " != " + format(rth)))
                # tup = (id, name, max_score, match_seq, rmax_score, rmatch_seq)
                match = MotifMatch(name, max_score, rmax_score, motif_index)
                matches.append(match)

        # output variant seq matches vs ref seq matches
        """for m in matches:
            line = "Variant matched motif "+m[1]+" ("+m[0]+") with score "
            line+= str(m[2])[:6]+" at seq: "+m[3]+"\nReference matched motif "+m[1]
            line+= "           with score "+str(m[4])[:6]+" at seq: "+m[5]
            #debug line+= "\n\tRefseq: "+ref_seq+"\n\tVarseq: "+var_seq
            print(line,file=output_f)"""    # WARNING: debug using undefined global

        return matches

    def motif_scores(self, baseline_p, sequence_str, wing_l, normalize):
        """
        Calculate if any motifs in the motif list match the given sequence.
        Requires that no motif have a length of 0.
        (former name: score_motifs, wanted grouped in class function list)

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')
            wing_l = length of sequence of bases flanking the variant
            normalize = boolean, if true then divide score by motif length
                concern: if not normalized, longer motifs can generate higher
                scores by length not true matches.

        Returns:
            matches = array of scores of the form:
                ( id,  name,  threshold,  max_score,  match_seq,  motif_index )

        see also: motif_scores_int
        """
        # list of matches between motifs and sequence
        scores = []

        for motif_index in range(len(self.motifs)):
            #### (iden, name, thresh, p_matrix) = tup
            motif_element = self.motifs[motif_index]

            # -- prepare the motif element and sequence string for scoring
            if not motif_element.valid_flag:
                continue

            if motif_element.matrix_type == 0:
                # need to calculate base occurence frequency by position. why
                # didn't it do this during insertion? will call many times here
                print("motif_scores: probabilities not previously calculated for " +
                    motif_element.print_str + ". Likely slower to calculate here. Fix")
                motif_element.calculate_probabilities()

            # trim flanking bases to one less than motif length
            # number of bases to trim
            trim_amount = wing_l - motif_element.positions + 1
            trim_seq = sequence.crop_from_left(sequence_str, trim_amount)
            trim_seq = sequence.crop_from_right(trim_seq, trim_amount)

            # -- actually do the scoring
            # highest match score for this motif (based on starting position)
            max_score = float("-inf")
            # sequence that matched motif best
            match_seq = ""

            # check match to all positions where motif overlaps with variant
            for pos in range(motif_element.positions):
                # check match starting at position
                # (score_motif will stop after the length of the motif)
                pos_score = motif_element.score_motif(baseline_p, trim_seq[pos:])
                if pos_score > max_score:
                    max_score = pos_score
                    match_seq = trim_seq[pos:pos + motif_element.positions]

            if normalize:
                max_score = max_score / motif_element.positions
                # YYY20161113.001: not sure if this is the correct way to handle threshold
                #    when normalizing the user may expect the value to be entered directly
                motif_threshold = motif_element.threshold / motif_element.positions
            else:
                motif_threshold = motif_element.threshold

            # debug print("Max match score:"+str(max_score)+" for motif "+name+" and
            # sequence "+match_seq+".")
            tupl = (motif_element.id, motif_element.name,
                    motif_threshold, max_score, match_seq, motif_index)
            scores.append(tupl)

        return scores

    def motif_scores_int(self, baseline_p, sequence_int, wing_l, normalize):
        """
        Calculate if any motifs in the motif list match the given sequence.
        Requires that no motif have a length of 0.
        (former name: score_motifs, wanted grouped in class function list)
        This version uses integer representations of the sequence

        #Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_int = array of characters (strings 'A','C','G', or 'T')
            wing_l = length of sequence of bases flanking the variant
            normalize = boolean, if true then divide score by motif length
                concern: if not normalized, longer motifs can generate higher
                scores by length not true matches.

        Returns:
            matches = array of scores of the form:
                ( id,  name,  threshold,  max_score,  match_seq,  motif_index )

        see also: motif_scores
        """
        # list of matches between motifs and sequence
        scores = []

        for motif_index in range(len(self.motifs)):
            #### (iden, name, thresh, p_matrix) = tup
            motif_element = self.motifs[motif_index]

            # -- prepare the motif element and sequence string for scoring
            if not motif_element.valid_flag:
                continue

            if motif_element.matrix_type == 0:
                # need to calculate base occurence frequency by position. why
                # didn't it do this during insertion? will call many times here
                print("motif_scores: probabilities not previously calculated for " +
                    motif_element.print_str + ". Likely slower to calculate here. Fix")
                motif_element.calculate_probabilities()

            # trim flanking bases to one less than motif length
            # number of bases to trim
            trim_amount = wing_l - motif_element.positions + 1
            trim_seq = sequence.crop_from_left(sequence_int, trim_amount)
            trim_seq = sequence.crop_from_right(trim_seq, trim_amount)

            # -- actually do the scoring
            # highest match score for this motif (based on starting position)
            max_score = float("-inf")
            # sequence that matched motif best
            match_seq = ""

            # check match to all positions where motif overlaps with variant
            # XXX: convolution has to be faster; then just pull peaks
            for pos in range(motif_element.positions):
                # check match starting at position
                # (score_motif will stop after the length of the motif)
                current_seq = trim_seq[pos:pos + motif_element.positions]
                pos_score = motif_element.score_motif_int(baseline_p, current_seq)
                if pos_score > max_score:
                    max_score = pos_score
                    match_seq = current_seq

            if normalize:
                max_score = max_score / motif_element.positions
                # YYY20161113.001: not sure if this is the correct way to handle threshold
                #    when normalizing the user may expect the value to be entered directly
                motif_threshold = motif_element.threshold / motif_element.positions
            else:
                motif_threshold = motif_element.threshold

            # debug print("Max match score:"+str(max_score)+" for motif "+name+" and
            # sequence "+match_seq+".")
            tupl = (motif_element.id, motif_element.name,
                    motif_threshold, max_score, match_seq, motif_index)
            scores.append(tupl)

        return scores

    def motif_scores_np(self, baseline_p, sequence_int, wing_l, normalize):
        """
        Calculate if any motifs in the motif list match the given sequence.
        Requires that no motif have a length of 0. ** This version uses numpy **
        (former name: score_motifs, wanted grouped in class function list)
        This version uses integer representations of the sequence passed as
        numpy arrays.

        #Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_int = array of characters (strings 'A','C','G', or 'T')
            wing_l = length of sequence of bases flanking the variant
            normalize = boolean, if true then divide score by motif length
                concern: if not normalized, longer motifs can generate higher
                scores by length not true matches.

        Returns:
            matches = array of scores of the form:
                ( id,  name,  threshold,  max_score,  match_seq,  motif_index )

        """
        # use baseline probability as np array with 0 element probability
        bp5 = np.insert(np.array(baseline_p), 0, 0)

        # list of matches between motifs and sequence
        scores = []

        for motif_index in range(len(self.motifs)):
            #### (iden, name, thresh, p_matrix) = tup
            motif_element = self.motifs[motif_index]

            # -- prepare the motif element and sequence string for scoring
            if not motif_element.valid_flag:
                continue

            if motif_element.matrix_type == 0:
                # need to calculate base occurence frequency by position. why
                # didn't it do this during insertion? will call many times here
                print("motif_scores: probabilities not previously calculated for " +
                    motif_element.print_str + ". Likely slower to calculate here. Fix")
                motif_element.calculate_probabilities()

            # trim flanking bases to one less than motif length
            # number of bases to trim
            trim_amount = wing_l - motif_element.positions + 1
            trim_seq = sequence.crop_from_left(sequence_int, trim_amount)
            trim_seq = sequence.crop_from_right(trim_seq, trim_amount)

            # -- actually do the scoring
            # highest match score for this motif (based on starting position)
            max_score = float("-inf")
            # sequence that matched motif best
            match_seq = ""

            # XXX: if switch to numpy version move this calculation to calculate_probabilities
            #    or when insert the element that way don't have to du for each score (repeated)
            # XXX: WARNING: RUNNING THIS COMMAND CORRUPTS MOTIF ELEMENT VARIABLE SPACE!
            #    ie the int and original uses of motif_element WILL NOT function properly
            #    makes the motif matrix np array and inserts a row of 0 elements
            #    0 elements used for not ACGT lookup, see sequence.convert2int
            motif_element.matrix = np.array(motif_element.matrix)
            motif_element.matrix = np.vstack(
                (np.repeat(0, motif_element.matrix.shape[1]), motif_element.matrix))

            # check match to all positions where motif overlaps with variant
            # XXX: convolution has to be faster; then just pull peaks
            for pos in range(motif_element.positions):
                # check match starting at position
                # (score_motif will stop after the length of the motif)
                current_seq = trim_seq[pos:pos + motif_element.positions]
                pos_score = motif_element.score_motif_np(bp5, current_seq)
                if pos_score > max_score:
                    max_score = pos_score
                    match_seq = current_seq

            if normalize:
                max_score = max_score / motif_element.positions
                # YYY20161113.001: not sure if this is the correct way to handle threshold
                #    when normalizing the user may expect the value to be entered directly
                motif_threshold = motif_element.threshold / motif_element.positions
            else:
                motif_threshold = motif_element.threshold

            # debug print("Max match score:"+str(max_score)+" for motif "+name+" and
            # sequence "+match_seq+".")
            tupl = (motif_element.id, motif_element.name,
                    motif_threshold, max_score, match_seq, motif_index)
            scores.append(tupl)

        return scores

    def process_local_env(self, baseline_p, matches, seq_element, co_binders_dict, v_seq, r_seq, wing_l):
        """
        Finds GC content, weak homotypic matches, and motif matches for
            co-binding transcription factors (not currently implemented).
        Homotypic match found by scoring wings only no overlap with variant

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            matches = list of MotifMatch objects, generated by motif_match()
            seq_element = SequenceElement object v_seq and r_seq belong to
            co_binders = dictionary that lists co-binding transcription factors for
                a given transcription factor name [feature not enabled. Use None]
            v_seq = variant sequence array of characters (strings 'A','C','G', or 'T')
            r_seq = reference sequence array of characters (strings 'A','C','G', or 'T')
            wing_l = integer wing length to override individual motif length

        Returns: Updates matches the list of MotifMatch objects

        XXX: WARNING: ERROR: this function was confusing class and other
            elements of the match
        XXX: this function never properly updated class elements;
            it is not clear what its intended functionality is/was
            review call usage & code generated based on assumptions/logic leaps.

        see also: process_local_env_int
        """

        # Get GC content
        # ZZZ: faster to convert to upper then count? Yes, actually, it is:
        # python -m timeit -s "vseq = 'GACCtcagTc'" "vseq.count('g') + vseq.count('G')"
        # 1000000 loops, best of 3: 0.384 usec per loop
        # python -m timeit -s "vseq = 'GACCtcagTc'" "vseq=vseq.upper()" "vseq.count('G')"
        # 1000000 loops, best of 3: 0.241 usec per loop
        v_seq = v_seq.upper()
        gc = v_seq.count('G') + v_seq.count('C')
        var_gc = gc / len(v_seq)

        r_seq = r_seq.upper()
        gc = r_seq.count('G') + r_seq.count('C')
        ref_gc = gc / len(r_seq)

        # Get list of motifs that matched
        for matches_index in range(len(matches)):
            # matches[matches_index]    # XXX
            match_index = matches[matches_index].motif_array_index
            if match_index < 0 or match_index > self.length():
                print(("**Error** process_local_env unable to find motif " +
                    matches[matches_index].name + " in MotifArray."))
                continue

            # motif info
            #self.motifs[match_index].matrix
            #self.motifs[match_index].threshold
            match_motif = self.motifs[match_index]

            # Find homotypic variant and reference matches
            vh_matches = []
            rh_matches = []

            # YYY!: original code and comments confusing:
            #    if homotopic matches are not supposed to overlap then...
            #    why does code in motifs.py do wing size 1 less than max positions?
            #    why did the original code below compute variant and reference
            #        matches against wings when if they don't overlap then same?
            #
            # ZZZ: which of the following should be used!?! (answer no overlap version)
            #   Check each wing separately,
            #   homotypic matches should not overlap the variant
            #   because variant and reference wings match only compute once
            #   wing size = size of the motif matched?
            # CCC-JA: Yeah, this part of the code was/is a mess. He was rushing towards the end
            #     to try to take into account homotypic matches and local GC content. I
            #     couldn't really decipher what was going on and he had a lot of unused variables
            #     and such left in here. It also didn't actually work quite right if I remember
            #     correctly. Our default wing size was +/- 50 bp from the variant position to
            #     look for homotypic matches, and yes, they should NOT overlap the variant.
            #     I suppose they might, however, overlap each other, but that'll get confusing
            #     in the output pretty quickly. If they overlap, best to just take the strongest
            #     match out of the overlaps and report that, I think.
            #
            #     GC content is pretty straightforward, just looking at % of bases in wings that are
            #     G or C in variant and reference sequences. Genuine TF binding sites tend to have
            #     slightly higher than normal GC content (usually ~40% in non-coding regions).
            #     If it was say, 55-60% in the local area of the variant, it might lend a bit of
            #     credence towards it being a genuine binding site. Just another piece of info.
            #     In summary, the original code was iffy, at best.
            # CCC-WK: code as written makes no attempt to check for wings overlapping; XXX

            # -- no overlap version
            left_wing = sequence.sub_from_end(seq_element.seq_left_wing.seq, wing_l)
                #match_motif.positions)
            right_wing = sequence.sub_from_start(seq_element.seq_right_wing.seq, wing_l)
                #match_motif.positions)
            vh_matches = match_motif.ht_matches_in(baseline_p, left_wing)
            vh_matches += match_motif.ht_matches_in(baseline_p, right_wing)
            # same strings so variant and ref homotypic matches are equal
            rh_matches = vh_matches

            # Update the match information
            matches[matches_index].var_ht = vh_matches
            matches[matches_index].var_gc = var_gc
            matches[matches_index].ref_ht = rh_matches
            matches[matches_index].ref_gc = ref_gc

        # Finding co-binders currently not implemented

        return matches

    def process_local_env_int(self, baseline_p, matches, seq_element, co_binders_dict, v_seq, r_seq, wing_l):
        """
        Integer sequence based version
        Finds GC content, weak homotypic matches, and motif matches for
            co-binding transcription factors (not currently implemented).
        Homotypic match found by scoring wings only no overlap with variant

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            matches = list of MotifMatch objects, generated by motif_match()
            seq_element = SequenceElement object v_seq and r_seq belong to
            co_binders = dictionary that lists co-binding transcription factors for
                a given transcription factor name [feature not enabled. Use None]
            v_seq = variant sequence array of integers (0-5)
            r_seq = reference sequence array of integers (0-5)
            wing_l = integer wing length to override individual motif length

        Returns: Updates matches the list of MotifMatch objects

        XXX: WARNING: ERROR: this function was confusing class and other
            elements of the match
        XXX: this function never properly updated class elements;
            it is not clear what its intended functionality is/was
            review call usage & code generated based on assumptions/logic leaps.

        see also: process_local_env
        """

        # Get GC content: counting the numeric equivalent see sequence.convert2int
        #    marginally faster as two separate lines than 1 line with () / len()
        gc = v_seq.count(2) + v_seq.count(3)
        var_gc = gc / len(v_seq)

        r_seq = r_seq.upper()
        gc = r_seq.count(2) + r_seq.count(3)
        ref_gc = gc / len(r_seq)

        # Get list of motifs that matched
        for matches_index in range(len(matches)):
            # matches[matches_index]    # XXX
            match_index = matches[matches_index].motif_array_index
            if match_index < 0 or match_index > self.length():
                print(("**Error** process_local_env unable to find motif " +
                    matches[matches_index].name + " in MotifArray."))
                continue

            # motif info
            #self.motifs[match_index].matrix
            #self.motifs[match_index].threshold
            match_motif = self.motifs[match_index]

            # Find homotypic variant and reference matches
            vh_matches = []
            rh_matches = []

            # YYY!: original code and comments confusing:
            #    if homotopic matches are not supposed to overlap then...
            #    why does code in motifs.py do wing size 1 less than max positions?
            #    why did the original code below compute variant and reference
            #        matches against wings when if they don't overlap then same?
            #
            # ZZZ: which of the following should be used!?! (answer no overlap version)
            #   Check each wing separately,
            #   homotypic matches should not overlap the variant
            #   because variant and reference wings match only compute once
            #   wing size = size of the motif matched?
            # CCC-JA: Yeah, this part of the code was/is a mess. He was rushing towards the end
            #     to try to take into account homotypic matches and local GC content. I
            #     couldn't really decipher what was going on and he had a lot of unused variables
            #     and such left in here. It also didn't actually work quite right if I remember
            #     correctly. Our default wing size was +/- 50 bp from the variant position to
            #     look for homotypic matches, and yes, they should NOT overlap the variant.
            #     I suppose they might, however, overlap each other, but that'll get confusing
            #     in the output pretty quickly. If they overlap, best to just take the strongest
            #     match out of the overlaps and report that, I think.
            #
            #     GC content is pretty straightforward, just looking at % of bases in wings that are
            #     G or C in variant and reference sequences. Genuine TF binding sites tend to have
            #     slightly higher than normal GC content (usually ~40% in non-coding regions).
            #     If it was say, 55-60% in the local area of the variant, it might lend a bit of
            #     credence towards it being a genuine binding site. Just another piece of info.
            #     In summary, the original code was iffy, at best.
            # CCC-WK: code as written makes no attempt to check for wings overlapping; XXX

            # -- no overlap version
            left_wing = sequence.sub_from_end(seq_element.seq_left_wing.seq, wing_l)
                #match_motif.positions)
            right_wing = sequence.sub_from_start(seq_element.seq_right_wing.seq, wing_l)
                #match_motif.positions)
            vh_matches = match_motif.ht_matches_in_int(baseline_p, left_wing)
            vh_matches += match_motif.ht_matches_in_int(baseline_p, right_wing)
            # same strings so variant and ref homotypic matches are equal
            rh_matches = vh_matches

            # Update the match information
            matches[matches_index].var_ht = vh_matches
            matches[matches_index].var_gc = var_gc
            matches[matches_index].ref_ht = rh_matches
            matches[matches_index].ref_gc = ref_gc

        # Finding co-binders currently not implemented

        return matches

    def search_set_max_position(self):
        """search elements in array to find and set max position information"""

        # set values to defaults
        self.max_positions = 0
        self.max_positions_index = -1

        # loop to set right value
        #for element in self.motifs:
        for index in range(0, (len(self.motifs) - 1)):
            if (self.motifs[index].positions > self.max_positions):
                self.max_positions = self.motifs[index].positions
                self.max_positions_index = index

        return


class MotifElement:
    """
    This class creates a Motif object
    Each object defines Motif characteristics principally the
    position matrix, matrix size, and modifications to the matrix
    """

    def __init__(self):
        self.name = "undefined"
        self.id = "undefined"
        self.matrix = [[] for x in range(4)]    # placeholder
            # matrix columns are positions
            # matrix rows are for bases, order ACGT
        self.matrix_type = 0    # 0 = count (Default), 1 = probability
                                # set by calculate_probabilities()
        self.positions = 0      # number of positions for motif matrix
                                # set by check_valid()
        self.valid_flag = False  # set by check_valid()
            # QQQ: store both types of matrix? memory hit is likely small
            # CCC-JA: Could, but probably not necessary since the counts aren't used again, right?
        self.threshold = 0
        self.base_probability = [] * 4   # base pair probability
        self.pseudocount = 0    # count number to add before probability calc

    def calculate_probabilities(self):
        """
        Given a count motif matrix, return motif probabilities
            matrix columns are positions
            matrix rows are for bases, order ACGT    XXX|QQQ make sure notes match

        Each input matrix element represents the number of times that base
        appears at that position in that motif for base defined by row and
        position defined by column.
        Each occurrences count is divided by the total number of occurrences
        for that position to give the frequency of each base by position.
        The position columns sum to 1.

        Pseudocounts are added in to avoid having a zero value in matrices with
        small sample sizes; see MotifElement.pseudocount
        """

        # only process count matrices that are valid
        if (self.matrix_type != 0 or not(self.valid_flag)):
            return

        # for clarity (also note: self.positions = length any matrix row)
        a = self.matrix[0]
        c = self.matrix[1]
        g = self.matrix[2]
        t = self.matrix[3]

        # QQQ: change to use numpy methods --> faster?
        # CCC-JA: Might be. Python's timeit function might be of use for testing
        #     small performance cases like this.
        # CCC-WK: the numpy reference is an extension of what I plan to do with convert2int
        #	the use not the conversion would get faster.
        # initialize output matrix of same size as input matrix
        new_m = [[0 for y in range(self.positions)] for x in range(4)]

        # divide each position by the sum of its column
        for pos in range(self.positions):
            total = float(a[pos]) + float(c[pos]) + float(g[pos]) + float(t[pos])
            total += 4 * self.pseudocount
            for base in range(4):
                new_m[base][pos] = (float(self.matrix[base][pos]) +
                     self.pseudocount) / total

        self.matrix = new_m
        self.matrix_type = 1

        return

    def check_valid(self):
        """
            sets valid_flag and positions for the MotifElement
            if all matrix elements are equal size --> set valid_flag = True
            if any matrix row length different than base length --> = False
        """
        self.valid_flag = True
        for mIndex in range(len(self.matrix)):
            if (mIndex == 0):
                self.positions = len(self.matrix[mIndex])
            else:
                if (self.positions != len(self.matrix[mIndex])):
                    self.valid_flag = False
                    break

        return

    def clear(self):
        """
        Resets current object;
        YYY: WARNING! python seems to allow this
            to affect previously assigned objects ie:
            b = MotifElement()
            c = MotifArray()
            c.add_motif(b)
            b.clear()    # also seems to clear info in c!

        """
        self.__init__()
        return

    def ht_matches_in(self, baseline_p, sequence_str):
        """
        Compute homotypic matches (scores above threshold) in the given sequence

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')
        Returns:

        see also: MotifArray.motif_scores()
        see also: ht_matches_in_int
        """
        homotypic_matches = []
        for pos in range(len(sequence_str) - self.positions + 1):
            pos_score = self.score_motif(baseline_p, sequence_str[pos:pos + self.positions])
            # Return the homotypic match if it is above the threshold
            # QQQ should this code normalize? see also YYY20161113.001
            if pos_score > self.threshold:
                homotypic_matches.append(pos_score)
        return homotypic_matches

    def ht_matches_in_int(self, baseline_p, sequence_int):
        """
        Compute homotypic matches (scores above threshold) in the given sequence

        Args:
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_int = array of integers (0-5)
        Returns:

        see also: MotifArray.motif_scores()
        see also: ht_matches_in
        """
        homotypic_matches = []
        for pos in range(len(sequence_int) - self.positions + 1):
            pos_score = self.score_motif_int(baseline_p, sequence_int[pos:pos + self.positions])
            # Return the homotypic match if it is above the threshold
            # QQQ should this code normalize? see also YYY20161113.001
            if pos_score > self.threshold:
                homotypic_matches.append(pos_score)
        return homotypic_matches

    def print_str(self):
        """return a string of the class object state"""
        print_string = self.name
        print_string += (":id:" + self.id)
        print_string += (":valid:" + self.valid_flag)
        print_string += (":matrix_type:" + self.matrix_type)
        print_string += (":positions:" + format(self.positions))
        print_string += (":threshold:" + format(self.threshold))
        print_string += (":pseudoCount:" + format(self.pseudocount))
        print_string += (":base_probability:" + format(self.base_probability))
        print_string += (":matrix:" + self.matrix)
        return print_string

    def score_base(self, probability, baseline_p):
        """
        compute the ... log likelihood ratio using log2

        Args:
            2 number inputs: probability, baseline_p
        Returns:
            log likelihood ration ( probability/baseline_p)

        YYY: note this could be outside the class  QQQ: why?
        """
        if (probability == 0):
            # should mathematically be negative infinity
            # this should not occur with pseudocounts added in
            return -100
        else:
            return log2(probability / baseline_p)

    def score_motif(self, baseline_p, sequence_str):
        """
        Calculate match between probability matrix for MotifElement and the
            given sequence_str starting at the beginning of the sequence
            and ending after the length of the motif.

        Args:
            p_matrix = probability matrix where rows are bases (A,C,G,T) and columns
                are probabilities (that sum to 1) at a given position. All rows should
                be the same length.
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_str = array of characters (strings 'A','C','G', or 'T')

        Returns:
            Match score between probability matrix and sequence_str

        see also: score_motif_int
        """

        score = 0

        if len(sequence_str) < self.positions:
            print("Sequence shorter than probability matrix. Returning score of 0.")
            return score

        if not self.valid_flag or self.matrix_type != 1:
            return score

        for pos in range(self.positions):
            # Match base: compute score for each overlap position
            if (sequence_str[pos] == 'A' or sequence_str[pos] == 'a'):
                score += self.score_base(self.matrix[0][pos], baseline_p[0])
            elif (sequence_str[pos] == 'C' or sequence_str[pos] == 'c'):
                score += self.score_base(self.matrix[1][pos], baseline_p[1])
            elif (sequence_str[pos] == 'G' or sequence_str[pos] == 'g'):
                score += self.score_base(self.matrix[2][pos], baseline_p[2])
            elif (sequence_str[pos] == 'T' or sequence_str[pos] == 't'):
                score += self.score_base(self.matrix[3][pos], baseline_p[3])
            # else is score = 0 therefore do nothing

        return score

    def score_motif_int(self, baseline_p, sequence_int):
        """
        Calculate match between probability matrix for MotifElement and the
            given sequence_int representation starting at the beginning of
            the sequence and ending after the length of the motif.

        Args:
            p_matrix = probability matrix where rows are bases (A,C,G,T) and columns
                are probabilities (that sum to 1) at a given position. All rows should
                be the same length.
            baseline_p = array of baseline probabilities of each base,
                in order (A, C, G, T). Probabilities should sum to 1.
                see get_baseline_probs()
            sequence_int = integers 0-5 for N ACGT

        Returns:
            Match score between probability matrix and sequence_str

        see also: score_motif
        """

        score = 0

        if len(sequence_int) < self.positions:
            print("Sequence shorter than probability matrix. Returning score of 0.")
            return score

        if not self.valid_flag or self.matrix_type != 1:
            return score

        for pos in range(self.positions):
            # Match base: compute score for each overlap position
            ind = sequence_int[pos] - 1
            if ind >= 0:
                score += self.score_base(self.matrix[ind][pos], baseline_p[ind])

        return score

    def score_motif_np(self, baseline_p, sequence_np):
        """
        Calculate match between probability matrix for MotifElement and the
            given sequence_int np array representation starting at the
            beginning of the sequence and ending after the length of the motif.
            This version uses integer representations of the sequence passed as
            numpy arrays.

        Args:
            p_matrix = numpy probability matrix where rows are bases (A,C,G,T) and columns
                are probabilities (that sum to 1) at a given position. All rows should
                be the same length. first row is all 0. p_matrix from self.matrix
            baseline_p = numpy array of baseline probabilities of each base,
                in order (N, A, C, G, T). Probabilities should sum to 1.
                N probability expected to be 0
                see get_baseline_probs()
            sequence_np = integers 0-5 for N ACGT

        Returns:
            Match score between probability matrix and sequence_str
        """

        score = 0

        if len(sequence_np) < self.positions:
            print("Sequence shorter than probability matrix. Returning score of 0.")
            return score

        if not self.valid_flag or self.matrix_type != 1:
            return score

        # Match base: compute cumulative score for each overlap position
        collookup = np.arange(self.positions)    # the position
            # collookup uses motif length not len(sequence_np)) should be equal
        xbp = baseline_p[sequence_np]            # base probability for ACGTN at position
        xvp = self.matrix[sequence_np, collookup]  # motif probability for ACGTN at position
        score = np.nansum((np.log2(xvp / xbp)))    # this replaces score_base

        return score


class MotifMatch:
    """
    This class could be rolled into MotifElement
    Current theory; faster not to, like normalized database table set,
    not duplicating MotifElement variables don't need and don't need to loop
    across the full MotifArray set
    """
    def __init__(self, name, var_score, ref_score, motif_index):
        # The transcription factor (tf) name
        self.name = name

        # tf motif match information (variant and reference sequence max log
        # likelihood score for the tf)
        self.var_score = var_score
        self.ref_score = ref_score

        # Does the motif match fall inside a chip peak for the same tf
        self.chip_match = False

        # Sequence environment data
        # variant and reference gc content
        self.var_gc = None
        self.ref_gc = None
        # variant and reference homotypic matches
        self.var_ht = []
        self.ref_ht = []    # XXX: drop see process_local_env notes
        # co-binding motifs (currently not implemented)
        self.cobinders = []

        self.motif_array_index = motif_index  # corresponding MotifArray object


# ______ START METHODS RELATED TO THE CLASS BUT NOT TIED TO IT _____
#    many could be made part of MotifArray or MotifElement

def get_motifs(motif_filename, pc, default_th, base_pr):
    """
    Read file containing the set of motif TF. Read from frequency matrix file.
    Computes and returns probability matrices for each motif TF except
    those individual motifs deemed to be invalid.

    Args:
        motif_f = filename of the frequency matrix
        pc = pseudocount value (small, greater than zero) to be added to all
            positions in matrix
        default_th = default threshold value used if none is listed in the
            motif file (should be third part of header if present).
            If None is given, default threshold will be calculated for motifs.
        base_pr = probabilities as a float array of the form:
            [ PrA, PrC, PrG, PrT ]
            ** Currently unused by this method
                but could be used to return a pssm instead of a pwm

    Returns:
        MotifArray object = set of Motif sequences as MotifElement objects
        Each object in the array defines Motif characteristics


        If reading multiple files use MotifArray join method on returned object

    YYY: The motif length for each TF is stored in MotifElement.positions
         The wing length of the compared sequence is by length of
         individual MotifElement not the max for all in the MotifArray

    """

    motif_set = MotifArray()
    bad_motif_count = 0    # error if count exceeds 3
    added_count = 0        # count of motifs added

    # Open and import motif file: note: with always closes open file
    with open(motif_filename) as file_handle:

        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n
        # arrays for C, G, T - each with same format as A
        # 1 blank line separates each individual TF position weight matrix
        # reading files assumes {header,A,C,G,T,blankline} order

        # index for line you are looking at in the motif matrix set (not file)
        i = 0

        # Iterate through file line by line.
        for line in file_handle:

            if line.startswith("#"):
                # QQQ: should the file read in a header # marked
                #      that specifies count or probability type?
                continue

            # First line contains id and name
            if i == 0:
                # must first remove brackets and split on whitespace
                line = line.strip().strip('[').strip(']').split()

                base_element = MotifElement()    # move creation to here
                                                 # prevents same object being tied to all
                                                 # prevents clear from affecting motifArray
                base_element.id = line[0]
                base_element.name = line[1]
                    #print(format(i) + ":" + base_element.id + " " + base_element.name)    # DEBUG LINE
                # If threshold given, use it for this matrix. Else, use default.
                if (len(line) > 2):
                    base_element.threshold = float(line[2])
                else:
                    base_element.threshold = default_th

            # Order of position weight matrices is A,C,G,T. Remove brackets etc.
            elif i < 5:
                # must first remove brackets and split on whitespace
                #     technically should modify index by found ACGT
                    #print(i)    # DEBUG LINE
                line = line.strip().strip('A').strip('C').strip('G').strip('T').strip().strip('[').strip(']').split()
                base_element.matrix[i - 1] = line[0:]

            # add motif to list and continue (there are 2 newlines between motifs)
            else:
                # check validity and set values passed
                base_element.check_valid()
                base_element.base_probability = base_pr
                base_element.pseudocount = pc

                if base_element.valid_flag:
                    # calculate base occurrence frequency by position
                    base_element.calculate_probabilities()
                else:
                    # invalid MotifElement matrix --> flag + do not include
                    # YYY: if invalid --> do not include;
                    #      file is suspect; if too many errors encountered
                    #      then throw error. may be safer to throw error if
                    #      any errors are encountered (ie bad file format)
                    bad_motif_count = bad_motif_count + 1
                    print(('Reading ' + format(motif_filename) +
                            ':\n\tdropping motif element ' +
                            format(base_element.name) + '.'))
                    if bad_motif_count > 2:
                        # XXX: throw an ERROR
                        print('Too many errors encountered. User input invalid')
                        raise ValueError('User defined motif file format invalid.')

                # append the current TF matrix to the motif set
                motif_set.add_motif(base_element)
                added_count = added_count + 1

                # reset the tracker variables
                i = -1
                #print(("before clear " + format(added_count) + ":" +
                #    format(motif_set.element_positions_list(True))))
            i += 1

    #print(("get_motifs read " + format(added_count) + " motifs from " + format(motif_filename)))
    #print(("before return set " + format(motif_set.element_positions_list(True))))
    return motif_set

def get_baseline_probs(baseline_f):
    """
    Read in baseline probabilities from a file name.

    Args:
        baseline_f a file containing a probability array of the form:
            [ PrA PrC PrG PrT ]
        Where PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)
        note1: the file can contain header lines marked by #
        note2: file format can technically be separated by arbitrary strings of
            any whitespace characters (space, tab, newline, return, formfeed).
            Only space tested during code development.
        note3: Separator must be whitespace, commas converted to empty string.
        note4: subsequent lines ignored once a properly formatted line is found

    Returns:
        Array with probabilities as a float array of the form:
        [ PrA, PrC, PrG, PrT ]
    """

    # Default baseline probability numbers (assumes all are equally likely)
    bp_array = [0.25, 0.25, 0.25, 0.25]

    # QQQ|YYY: note safer if check for invalid files everywhere
    # YYY: Yeah, a function to just validate the input file formats is likely a good idea.
    with open(baseline_f) as f:
        try:
            for line in f:
                # remove header
                if line.startswith("#"):
                    continue
                # remove commas, brackets, and whitespace on far left and right
                line = line.strip().replace('[', '').replace(']', '').replace(',', '')
                if line != "":
                    line = line.split()
                    if (len(line) >= 4):
                        for idx in range(4):
                            bp_array[idx] = float(line[idx])
                        return bp_array
        except ValueError:
            print(("**ERROR** Baseline probability file incorrectly formatted.\n" +
                  "\tFile should contain only [ PrA PrC PrG PrT ] \n" +
                  "\tWhere PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)\n" +
                  "\tContinuing with default probabilities: " + format(bp_array)))
            return bp_array

        print(("**ERROR** Empty file found.\n" +
              "\tFile should contain only [ PrA PrC PrG PrT ] \n" +
              "\tWhere PrA + PrC + PrG + PrT = 1 (and all are positive and non-zero)\n" +
              "\tContinuing with default probabilities: " + format(bp_array)))
        return bp_array


def log2infforced(value):
    """
    compute the log2 with <= 0 forced to -100

    Args:
        value
    Returns:
        log2 or -100 if <= 0 input
    """

    if (value <= 0):
        # avoids mathematically negative infinity
        return -100
    else:
        return log2(value)


def log2infforcednp(nparray):
    """
    call log2infforced for first element in the array

    Args:
        array
    Returns:
        log2infforced(nparray([0]))
    """

    return log2infforced(nparray[0])


def put_motifs(motif_filename, motif_set, pc, default_th, base_pr):
    """
    Write file containing the set of motif TF. Read from frequency matrix file.
    Computes and returns probability matrices for each motif TF except
    those individual motifs deemed to be invalid.

    XXX: function is incomplete

    Args:
        motif_f = filename of the frequency matrix
        pc = pseudocount value (small, greater than zero) to be added to all
            positions in matrix
        default_th = default threshold value used if none is listed in the
            motif file (should be third part of header if present).
            If None is given, default threshold will be calculated for motifs.
        base_pr = probabilities as a float array of the form:
            [ PrA, PrC, PrG, PrT ]
            ** Currently unused by this method
                but could be used to return a pssm instead of a pwm

    Returns:
        MotifArray object = set of Motif sequences as MotifElement objects
        Each object in the array defines Motif characteristics


        If reading multiple files use MotifArray join method on returned object

    YYY: The motif length for each TF is stored in MotifElement.positions
         The wing length of the compared sequence is by length of
         individual MotifElement not the max for all in the MotifArray

    """
    return -1
    motif_set = MotifArray()
    bad_motif_count = 0     # error if count exceeds 3

    # Open and import motif file: note: with always closes open file
    with open(motif_filename) as file_handle:

        # JASPAR motif file has >name \n A [ tab delineated weight array ] \n
        # arrays for C, G, T - each with same format as A
        # 1 blank line separates each individual TF position weight matrix
        # reading files assumes {header,A,C,G,T,blankline} order

        # index for line you are looking at in the motif matrix set (not file)
        i = 0

        # Iterate through file line by line.
        for line in file_handle:

            # First line contains id and name
            if i == 0:
                # must first remove brackets and split on whitespace
                line = line.strip().strip('[').strip(']').split()

                base_element = MotifElement()
                base_element.id = line[0]
                base_element.name = line[1]
                    #print(format(i) + ":" + base_element.id + " " + base_element.name)    # DEBUG LINE
                # If threshold given, use it for this matrix. Else, use default.
                if (len(line) > 2):
                    base_element.threshold = float(line[2])
                else:
                    base_element.threshold = default_th

            # Order of position weight matrices is A,C,G,T. Remove brackets etc.
            elif i < 5:
                # must first remove brackets and split on whitespace
                #     technically should modify index by found ACGT
                    #print(i)    # DEBUG LINE
                line = line.strip().strip('A').strip('C').strip('G').strip('T').strip().strip('[').strip(']').split()
                base_element.matrix[i - 1] = line[0:]

            # add motif to list and continue (there are 2 newlines between motifs)
            else:
                # check validity and set values passed
                base_element.check_valid()
                base_element.base_probability = base_pr
                base_element.pseudocount = pc

                if base_element.valid_flag:
                    # calculate base occurrence frequency by position
                    base_element.calculate_probabilities()
                else:
                    # invalid MotifElement matrix --> flag + do not include
                    # YYY: if invalid --> do not include;
                    #      file is suspect; if too many errors encountered
                    #      then throw error. may be safer to throw error if
                    #      any errors are encountered (ie bad file format)
                    bad_motif_count = bad_motif_count + 1
                    print('Reading ' + format(motif_filename) +
                            ':\n\tdropping motif element ' +
                            format(base_element.name) + '.' )
                    if bad_motif_count > 2:
                        # XXX: throw an ERROR
                        print('Too many errors encountered. User input invalid')
                        raise ValueError('User defined motif file format invalid.')

                # append the current TF matrix to the motif set
                motif_set.add_motif(base_element)

                # reset the tracker variables
                i = -1
                base_element.clear()    # BAD IDEA
            i += 1

    return motif_set
