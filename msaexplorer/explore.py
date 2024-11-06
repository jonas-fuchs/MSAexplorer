"""
This contains ...
"""

# built-in
import math
import re
# installed
from Bio import AlignIO
import numpy as np
from numpy import ndarray
# msaexplorer
from msaexplorer import config


class Alignment:
    """
    An alignment class that allows computation of several stats
    """

    def __init__(self, alignment_path: str, reference_id:str=None, zoom_range:tuple=None):
        """
        Initialise an Alignment object.
        :param alignment_path: path to alignment file
        :param reference_id: reference id
        :param zoom_range: start and stop positions to zoom into the alignment
        """
        self._alignment = Alignment._read_alignment(alignment_path)
        self._reference_id = self._validate_ref(reference_id, self._alignment)
        self._zoom = self._validate_zoom(zoom_range, self._alignment)
        self._aln_type = self._determine_aln_type(self._alignment)

    # Static methods
    @staticmethod
    def _read_alignment(alignment_path: str) -> dict:
        """
        Read alignment with Bio.AlignIO.
        :param alignment_path: path to alignment file
        :return: dictionary with ids as keys and sequences as values
        """
        alignment = dict()

        for sequence in AlignIO.read(alignment_path, "fasta"):
            alignment[sequence.id] = str(sequence.seq)

        return alignment

    @staticmethod
    def _validate_ref(reference: str|None, alignment: dict) -> str|None|ValueError:
        """
        Validate if the ref seq is indeed part of the alignment.
        :param reference: reference seq id
        :param alignment: alignment dict
        :return: validated reference
        """
        if reference in alignment.keys():
            return reference
        elif reference is None:
            return reference
        else:
            raise ValueError('Reference not in alignment.')

    @staticmethod
    def _validate_zoom(zoom: tuple, original_aln:dict) -> ValueError|tuple|None:
        """
        Validates if the user defined zoom range is within the start, end of the initial
        alignment.\n
        :param zoom: zoom range
        :param original_aln: non-zoomed alignment dict
        :return: validated zoom range
        """
        if zoom is not None:
            if len(zoom) != 2:
                raise ValueError('Zoom position have to be (zoom_start, zoom_end)')
            for position in zoom:
                if type(position) != int:
                    raise ValueError('Zoom positions have to be integers.')
                if position not in range(0, len(original_aln[list(original_aln.keys())[0]])):
                    raise ValueError('Zoom position out of range')

        return zoom

    @staticmethod
    def _shannons_entropy(chars: list, states: int, aln_type: str) -> float:
        """
        Calculate the shannon's entropy of a sequence and
        normalized between 0 and 1.
        :param chars: characters at an alignment position
        :param states: number of potential characters that can be present
        :returns: entropy
        """
        ent, n_chars = 0, len(chars)
        # only one char is in the list
        if n_chars <= 1:
            return ent
        # calculate the number of unique chars and their counts
        chars, char_counts = np.unique(chars, return_counts=True)
        char_counts = char_counts.astype(float)
        # ignore gaps for entropy calc
        char_counts, chars = char_counts[chars != "-"], chars[chars != "-"]
        # correctly handle ambiguous nucleotides
        if states == 4:
            index_to_drop = []
            for index, char in enumerate(chars):
                if char in config.AMBIG_NUCS[aln_type]:
                    index_to_drop.append(index)
                    nucleotide, nt_counts = np.unique(config.AMBIG_NUCS[aln_type][char], return_counts=True)
                    nt_counts = nt_counts / len(config.AMBIG_NUCS[aln_type][char])
                    # add the proportionate numbers to initial array
                    for nucleotide, nt_count in zip(nucleotide, nt_counts):
                        if nucleotide in chars:
                            char_counts[chars == nucleotide] += nt_count
                        else:
                            chars, char_counts = np.append(chars, nucleotide), np.append(char_counts, nt_count)
            # drop the ambiguous characters from array
            char_counts, chars = np.delete(char_counts, index_to_drop), np.delete(chars, index_to_drop)
        # correctly handle 'N' in as:
        if states == 20 and 'N' in chars:
            temp_counts =  char_counts[chars == 'N']/states
            char_counts, chars = char_counts[chars != 'N'], chars[chars != 'N']
            char_counts += temp_counts
            for amino_acid in config.amino_acids:
                if amino_acid in chars:
                    continue
                chars, char_counts = np.append(chars, amino_acid), np.append(char_counts, temp_counts)
        # calc the entropy
        probs = char_counts / n_chars
        if np.count_nonzero(probs) <= 1:
            return ent
        for prob in probs:
            ent -= prob * math.log(prob, states)

        return ent

    @staticmethod
    def _determine_aln_type(alignment) -> str:
        """
        Determine the most likely type of alignment
        if 70% of chars in the alignment are nucleotide
        chars it is most likely a nt alignment
        :return: type of alignment
        """
        counter = int()
        for record in alignment:
            if 'U' in record:
                return 'RNA'
            counter += sum(map(alignment[record].count, ['A', 'C', 'G', 'T', 'N', '-']))
        # determine which is the most likely type
        if counter/len(alignment) >= 0.7 * len(alignment[list(alignment.keys())[0]]):
            return 'DNA'
        else:
            return 'AS'

    # Properties with setters
    @property
    def reference_id(self):
        return self._reference_id

    @reference_id.setter
    def reference_id(self, ref_id: str):
        """
        Set and validate the reference id.
        """
        self._reference_id = self._validate_ref(ref_id, self.alignment)

    @property
    def zoom(self) -> tuple:
        return self._zoom

    @zoom.setter
    def zoom(self, zoom_pos: tuple):
        """
        Validate if the user defined zoom range.
        """
        self._zoom = self._validate_zoom(zoom_pos, self._alignment)

    # Property without setters
    @property
    def aln_type(self) -> str:
        """
        define the aln type:
        RNA, DNA or AS
        """
        return self._aln_type

    # On the fly properties without setters
    @property
    def length(self) -> int:
        """
        Determine the length of the (sliced) alignment.
        """
        for record in self.alignment:
            return len(self.alignment[record])

    @property
    def alignment(self) -> dict:
        """
        (zoomed) version of the alignment.
        """
        if self.zoom is not None:
            zoomed_aln = dict()
            for seq in self._alignment:
                zoomed_aln[seq] = self._alignment[seq][self.zoom[0]:self.zoom[1]]
            return zoomed_aln
        else:
            return self._alignment

    # functions for different alignment stats
    def get_reference_coords(self) -> tuple:
        """
        Determine the start and end coordinates of the reference sequence
        defined as the first/last nucleotide in the reference sequence
        (excluding N and gaps).

        :return: start, end
        """
        start, end = 0, self.length

        if self.reference_id is None:
            return start, end
        else:
            # 5' --> 3'
            for start in range(self.length):
                if self.alignment[self.reference_id][start] not in ['-','N']:
                    break
            # 3' --> 5'
            for end in range(self.length-1, 0, -1):
                if self.alignment[self.reference_id][end] not in ['-','N']:
                    break

            return start, end

    def calc_length_stats(self) -> dict:
        """
        Determine the stats for the length of the ungapped seqs in the alignment.
        :return:
        """

        seq_lengths = [len(self.alignment[x].replace('-', '')) for x in self.alignment]

        return {'number of seq': len(self.alignment),
                'mean length': float(np.mean(seq_lengths)),
                'std length': float(np.std(seq_lengths)),
                'min length': int(np.min(seq_lengths)),
                'max length': int(np.max(seq_lengths))
                }

    def calc_entropy(self) -> list:
        """
        Calculate the entropy for every position in an alignment.
        :return: Entropies at each position.
        """

        aln = self.alignment
        entropys = []

        if self.aln_type == 'AS':
            states = 20
        else:
            states = 4
        # iterate over alignment positions and the sequences
        for nuc_pos in range(self.length):
            pos = []
            for record in aln:
                pos.append(aln[record][nuc_pos])
            entropys.append(Alignment._shannons_entropy(pos, states, self.aln_type))

        return entropys

    def calc_gc(self) -> list|TypeError:
        """
        Determine the GC content for every position in an nt alignment.
        :return: GC content for every position.
        """
        if self.aln_type == 'AS':
            raise TypeError("GC computation is not possible for aminoacid alignment")

        gc, aln, amb_nucs = [], self.alignment, config.AMBIG_NUCS[self.aln_type]

        for position in range(self.length):
            nucleotides = str()
            for record in aln:
                nucleotides = nucleotides + aln[record][position]
            # ini dict with chars that occur and which ones to
            # count in which freq
            to_count = {
                'G': 1,
                'C': 1,
                '-': 0.5  # similar to 'N' its not clear what the freq is at a deletion
            }
            # handle ambig. nuc
            for char in amb_nucs:
                if char in nucleotides:
                    to_count[char] = (amb_nucs[char].count('C') + amb_nucs[char].count('G'))/len(amb_nucs[char])

            gc.append(
                sum([nucleotides.count(x)*to_count[x] for x in to_count])/len(nucleotides)
            )

        return gc

    def calc_coverage(self) -> list:
        """
        Determine the coverage of every position in an alignment.
        This is defined as 1 - cumulative length of '-' characters
        :return: Coverage at each alignment position.
        """
        coverage, aln = [], self.alignment

        for nuc_pos in range(self.length):
            pos = str()
            for record in aln:
                pos = pos + aln[record][nuc_pos]
            coverage.append(1-pos.count('-')/len(pos))

        return coverage
    # TODO more consensus options?
    def get_consensus(self) -> str:
        """
        Creates a non-gapped consensus sequence with the most freq
        nucleotide.
        :return: consensus sequence
        """
        consensus, aln = "", self.alignment

        for nuc_pos in range(self.length):
            chars = []
            for record in aln:
                chars.append(aln[record][nuc_pos])
            chars, values = np.unique(chars, return_counts=True)
            chars, values = chars[chars != '-'], values[chars != '-']
            possible_chars = chars[values == values.max()]
            # if max function yields two or more values ensure that the one
            # specified for the consensus is not an 'N'
            if len(possible_chars) > 1:
                possible_chars = possible_chars[possible_chars != 'N']
            consensus = consensus + str(possible_chars[0])

        return consensus

    def calc_reverse_complement_alignment(self) -> dict|TypeError:
        """
        Reverse complement the alignment.
        :return: Alignment
        """
        if self.aln_type == 'AS':
            raise TypeError('Reverse complement only for RNA or DNA.')

        aln = self.alignment
        reverse_complement_dict = {}

        for seq_id in aln:
            reverse_complement_dict[seq_id] = ''.join(config.complement[base] for base in reversed(aln[seq_id]))

        return reverse_complement_dict

    def calc_identity_alignment(self) -> ndarray:
        """
        Converts alignment to identity array compared to consensus or reference:\n
        nan: deletion \n
        0: identical \n
        1: mismatch
        :return: identity alignment
        """

        aln = self.alignment

        if self.reference_id is not None:
            ref = aln[self.reference_id]
        else:
            ref = self.get_consensus()

        identity_matrix, aln = [], self.alignment

        for seq_id in aln:
            temp_list = []
            for i, char in enumerate(aln[seq_id]):
                if char == '-':
                    temp_list.append(np.nan)
                elif char == ref[i]:
                    temp_list.append(0)
                else:
                    temp_list.append(1)
            identity_matrix.append(temp_list)

        return np.array(identity_matrix, dtype=object)
    # TODO write the function
    def calc_similarity_alignment(self):
        """
        Calculate the similarity score between the alignment and the reference seq based on
        Similarity matrixes (Blossum etc)
        :return:
        """
        pass

    def calc_percent_recovery(self) -> dict:
        """
        Recovery per sequence either compared to a consensus seq
        or the reference seq.\n
        Defined as:\n
        100 - (number of non-N and non-gap characters of non-gapped reference regions) / length of ungapped region * 100
        :return: dict
        """

        aln = self.alignment

        if self.reference_id is not None:
            ref = aln[self.reference_id]
        else:
            ref = self.get_consensus()

        # define positions in reference that are not gapped
        non_gaps = [(match.start(), match.end()) for match in re.finditer(r"[^-]+", ref)]
        cumulative_length = sum(end - start for start, end in non_gaps)

        # count 'N' and '-' chars in non-gapped regions
        recovery_over_ref = dict()

        for seq_id in aln:
            if seq_id == self.reference_id:
               continue
            recovery_over_ref[seq_id] = 0
            for region in non_gaps:
                recovery_over_ref[seq_id] += aln[seq_id][region[0]:region[1]].count('-')
                recovery_over_ref[seq_id] += aln[seq_id][region[0]:region[1]].count('N')
            recovery_over_ref[seq_id] = 100 - recovery_over_ref[seq_id] / cumulative_length * 100

        return recovery_over_ref

    def calc_character_frequencies(self) -> dict:
        """
        Calculate the percentage characters in the alignment:
        The frequencies are counter by seq and in total. The
        percentage of non-gap characters in the alignment is
        relative to the total number of non-gap characters.
        The gap percentage is relative to the sequence length.

        The output is a nested dictionary.

        :return: Character frequencies
        """

        aln, aln_length = self.alignment, self.length

        freqs = {'total':{'-': {'counts':0, '% of alignment':float()}}}

        for seq_id in aln:
            freqs[seq_id], all_chars = {'-': {'counts':0, '% of alignment':float()}}, 0
            unique_chars = set(aln[seq_id])
            for char in unique_chars:
                if char == '-':
                    continue
                # add characters to dictionaries
                if char not in freqs[seq_id]:
                    freqs[seq_id][char] = {'counts':0, '% of non-gapped':0}
                if char not in freqs['total']:
                    freqs['total'][char] = {'counts':0, '% of non-gapped':0}
                # count non-gap chars
                freqs[seq_id][char]['counts'] += aln[seq_id].count(char)
                freqs['total'][char]['counts'] += freqs[seq_id][char]['counts']
                all_chars += freqs[seq_id][char]['counts']
            # normalize counts
            for char in freqs[seq_id]:
                if char == '-':
                    continue
                freqs[seq_id][char]['% of non-gapped'] = freqs[seq_id][char]['counts'] / all_chars * 100
                freqs['total'][char]['% of non-gapped'] += freqs[seq_id][char]['% of non-gapped']
            # count gaps
            freqs[seq_id]['-']['counts'] = aln[seq_id].count('-')
            freqs['total']['-']['counts'] += freqs[seq_id]['-']['counts']
            # normalize gap counts
            freqs[seq_id]['-']['% of alignment'] = freqs[seq_id]['-']['counts'] / aln_length * 100
            freqs['total']['-']['% of alignment'] += freqs[seq_id]['-']['% of alignment']

        # normalize the total counts
        for char in freqs['total']:
            for value in freqs['total'][char]:
                if value == '% of alignment' or value == '% of non-gapped':
                    freqs['total'][char][value] = freqs['total'][char][value] / len(aln)

        return freqs
    # TODO write the function
    def calc_pairwise_distance(self):
        """
        score for each seq how many identical characters it has compared to each other seq
        :return:
        """
        pass
    # TODO fix positions for reverse complement
    def get_conserved_orfs(self, min_length:int=500) -> dict:
        """
        conserved ORF definition:
            - conserved starts and stops
            - start, stop must be on the same frame
            - stop - start must be at least min_length
            - all ungapped seqs[start:stop] must have at least min_length
            - no ungapped seq can have a Stop in between Start Stop

        Algorithm overview:
            - check for conserved start and stop codons
            - iterate over all three frames
            - check each start and next sufficiently far away stop codon
            - check if all ungapped seqs between start and stop codon are >= min_length
            - check if no ungapped seq in the alignment has a stop codon
            - write to dictionary
            - classify as internal if the stop codon has already been written with a prior start
            - repeat for reverse complement

        :return: ORF positions and internal ORF positions
        """

        if self.aln_type == 'AS':
            raise TypeError('ORF search only for RNA/DNA alignments')

        starts = config.start_codons[self.aln_type]
        stops = config.stop_codons[self.aln_type]
        seq = self.alignment[list(self.alignment.keys())[0]]
        identities = self.calc_identity_alignment()
        alignments = [self.alignment, self.calc_reverse_complement_alignment()]
        aln_len = self.length

        # ini
        orf_counter = 0
        orf_dict = {}

        for aln, direction in zip(alignments, ['+', '-']):
            # check for starts and stops in the first seq and then check if these are present in all seqs
            conserved_starts, conserved_stops = [], []
            for pos in range(aln_len):
                if seq[pos:pos+3] in starts and np.nansum(identities[:,[x for x in range(pos,pos+3)]]) == 0:
                    conserved_starts.append(pos)
                if seq[pos:pos+3] in stops and np.nansum(identities[:,[x for x in range(pos,pos+3)]]) == 0:
                    conserved_stops.append(pos)
            # check each frame
            for frame in (0, 1, 2):
                potential_starts = [x for x in conserved_starts if x % 3 == frame]
                potential_stops = [x for x in conserved_stops if x % 3 == frame]
                last_stop = -1
                for start in potential_starts:
                    # go to the next stop that is sufficiently far away in the alignment
                    next_stops = [x for x in potential_stops if x+2 >= start + min_length]
                    if not next_stops:
                        continue
                    next_stop = next_stops[0]
                    # re-check the lengths of all ungapped seqs
                    ungapped_seq_lengths = [len(aln[x][start:next_stop+2].replace('-', '')) >= min_length for x in aln.keys()]
                    if not all(ungapped_seq_lengths):
                        continue
                    # set checker if any seq has an additional stop
                    additional_stop = False
                    # re-check every single ungapped seq for the presence of an in-frame stop
                    for seq_id in aln:
                        sliced_seq = aln[seq_id][start:next_stop].replace('-', '')
                        for pos in range(0, len(sliced_seq), 3):
                            if sliced_seq[pos:pos + 3] in stops:
                                additional_stop = True
                                break
                        if additional_stop:
                            break
                    # if no stop codon between start and stop --> write to dictionary
                    if not additional_stop:
                        if last_stop != next_stop:
                            if direction == '+':
                                positions = [start, next_stop + 3]
                            else:
                                positions = [aln_len - next_stop - 3, aln_len - start]

                            orf_dict[f'ORF_{orf_counter}'] = {'positions': positions,
                                                              'frame': frame,
                                                              'strand': direction,
                                                              'internal': []
                                                              }
                            orf_counter += 1
                            last_stop = next_stop
                        else:
                            orf_dict[f'ORF_{orf_counter - 1}']['internal'].append((start, next_stop))

        return orf_dict


def read_annotation():
    pass


def parse_bed():
    pass


def parse_gff():
    pass


def parse_genbank():
    pass


# functions to save data in standard formats
def dict_to_bed_file():
    pass