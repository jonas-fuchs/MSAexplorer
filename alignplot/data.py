import math
# import installed
from Bio import AlignIO
import numpy as np
# alignplot
from alignplot import config


class Alignment:
    """
    an alignment class that allows computation of several stats
    """

    def __init__(self, init_alignment, reference=None):
        self.alignment = Alignment.read_alignment(init_alignment)
        self._reference = reference

    @staticmethod
    def shannons_entropy(chars, states):
        """
        input is a list of as or nt characters
        """
        ent, n_chars = 0, len(chars)
        # only one char is in the list
        if n_chars <= 1:
            return ent
        # calculate the number of unique chars and their counts
        values, counts = np.unique(chars, return_counts=True)
        counts = counts.astype(float)
        # ignore gaps
        counts, values = counts[values != "-"], values[values != "-"]
        # correctly handle ambiguous nucleotides
        if states == 4:
            to_drop = []
            for i, value in enumerate(values):
                if value in config.AMBIG_NUCS:
                    to_drop.append(i)
                    amb_values, amb_counts = np.unique(config.AMBIG_NUCS[value], return_counts=True)
                    amb_counts = amb_counts / len(config.AMBIG_NUCS[value])
                    # add the proportionate numbers to initial array
                    for amb_value, amb_count in zip(amb_values, amb_counts):
                        if amb_value in values:
                            counts[values == amb_value] += amb_count
                        else:
                            values, counts = np.append(values, amb_value), np.append(counts, amb_count)
            # drop the ambiguous characters from array
            counts, values = np.delete(counts, to_drop), np.delete(values, to_drop)
        # correctly handle 'N' in as:
        if states == 20 and 'N' in values:
            temp_counts =  counts[values == 'N']/states
            counts, values = counts[values != 'N'], values[values != 'N']
            counts += temp_counts
            for amino_acid in config.amino_acids:
                if amino_acid in values:
                    continue
                values, counts = np.append(values, amino_acid), np.append(counts, temp_counts)
        # calc the entropy
        probs = counts / n_chars
        if np.count_nonzero(probs) <= 1:
            return ent
        for prob in probs:
            ent -= prob * math.log(prob, states)

        return ent

    @staticmethod
    def read_alignment(alignment_path):
        """
        read alignment with AlignIO and
        convert to dict
        """
        alignment = dict()

        for sequence in AlignIO.read(alignment_path, "fasta"):
            alignment[sequence.id] = str(sequence.seq)

        return alignment

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, value):
        if value in self.alignment.keys() or value is None:
            self._reference = value
        else:
            raise ValueError('Reference not in alignment and not None')

    @property
    def type(self):
        """
        determine the most likely type of alignment
        :return: type of alignment
        """
        counter = int()
        for record in self.alignment:
            counter += sum(map(self.alignment[record].count, ['A', 'C', 'G', 'T']))
        # determine which is the most likely type
        if counter/len(self.alignment) > 0.7:
            return 'nt'
        else:
            return 'as'

    @property
    def length(self):
        """
        determine the length of the alignment
        :return: length of alignment
        """
        for record in self.alignment:
            return len(self.alignment[record])

    @property
    def entropy(self):
        """
        calculate the entropy for every position in an alignment.
        :return: entropy list
        """
        entropys = []
        if self.type == 'as':
            states = 20
        else:
            states = 4
        # iterate over alignment positions and the sequences
        for nuc_pos in range(0, self.length):
            pos = []
            for record in self.alignment:
                pos.append(self.alignment[record][nuc_pos])
            entropys.append(Alignment.shannons_entropy(pos, states))

        return entropys

    @property
    def gc(self):
        """
        determine the GC content for every position in an alignment.
        :return: GC list
        """
        if self.type == 'as':
            raise TypeError("gc computation is not possible for aminoacid alignment")

        gc = []
        for nuc_pos in range(0, self.length):
            pos = str()
            for record in self.alignment:
                pos = pos + self.alignment[record][nuc_pos]
            # ini dict with chars that occur and which ones to
            # count in which freq
            to_count = {
                'G': 1,
                'C': 1,
                '-': 0.5  # similar to 'N' its not clear what the freq is at a deletion
            }
            # handle ambig nuc
            for char in config.AMBIG_NUCS:
                if char in pos:
                    to_count[char] = (config.AMBIG_NUCS[char].count('C') + config.AMBIG_NUCS[char].count('G'))/len(config.AMBIG_NUCS[char])

            gc.append(
                sum([pos.count(x)*to_count[x] for x in to_count])/len(pos)
            )

        return gc

    @property
    def coverage(self):
        """
        determine the coverage of every position in an alignment.
        defined as 1 - cumulative length of '-' characters
        :return: coverage list
        """
        coverage = []
        for nuc_pos in range(0, self.length):
            pos = str()
            for record in self.alignment:
                pos = pos + self.alignment[record][nuc_pos]
            coverage.append(1-pos.count('-')/len(pos))

        return coverage

    @property
    def consensus(self):
        """
        creates a consensus sequence with the most freq
        nucleotide
        :return: consensus as string
        """
        consensus = ""
        for nuc_pos in range(0, self.length):
            chars = []
            for record in self.alignment:
                chars.append(self.alignment[record][nuc_pos])
            chars, values = np.unique(chars, return_counts=True)
            chars, values = chars[chars != '-'], values[chars != '-']
            possible_chars = chars[values == values.max()]
            # if max function yields two or more values ensure that the one
            # specified for the consensus is not an 'N'
            if len(possible_chars) > 1:
                possible_chars = possible_chars[possible_chars != 'N']
            consensus = consensus + str(possible_chars[0])

        return consensus

    @property
    def identity_array(self):
        """
        converts alignment to identity array
        either compared to a consensus seq
        harboring the most frequent nucleotide
        or a sepcified reference seq
        'N' = -2
        gap = -1
        identical = 0
        mismatch = 1
        :return: identity array
        """
        if self.reference is not None:
            ref = self.alignment[self.reference]
        else:
            ref = self.consensus

        identity_matrix = []
        for seq_id in self.alignment:
            if seq_id == self.reference:
                continue
            temp_list = []
            for i, nuc in enumerate(self.alignment[seq_id]):
                if nuc == '-':
                    temp_list.append(-1)
                elif nuc == 'N' and ref[i] != 'N':
                    temp_list.append(-2)
                elif nuc == ref[i]:
                    temp_list.append(0)
                else:
                    temp_list.append(1)
            identity_matrix.append(temp_list)

        return np.array(identity_matrix)

    #@property
    def mean_recovery(self):
        """
        mean recovery per sequence either compared to a consensus seq
        or the reference seq. defined as:
        (alignment_length - (total gaps + number of amb. char))/(ref_length)
        ????
        :return:
        """
        pass