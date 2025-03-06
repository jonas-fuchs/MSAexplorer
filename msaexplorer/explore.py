"""
# Explore module

This module contains the two classes `MSA` and `Annotation` which are used to read in the respective files
and can be used to compute several statistics or be used as the input for the `draw` module functions.

## Classes

"""

# built-in
import math
import collections
import re
from copy import deepcopy
from typing import Callable, Dict

# installed
import numpy as np
from numpy import ndarray

# msaexplorer
from msaexplorer import config

# TODO: test with as alignment
class MSA:
    """
    An alignment class that allows computation of several stats
    """

    def __init__(self, alignment_path: str, reference_id: str = None, zoom_range: tuple = None):
        """
        Initialise an Alignment object.
        :param alignment_path: path to alignment file
        :param reference_id: reference id
        :param zoom_range: start and stop positions to zoom into the alignment
        """
        self._alignment = self._read_alignment(alignment_path)
        self._reference_id = self._validate_ref(reference_id, self._alignment)
        self._zoom = self._validate_zoom(zoom_range, self._alignment)
        self._aln_type = self._determine_aln_type(self._alignment)

    # Static methods
    @staticmethod
    def _read_alignment(file_path: str) -> dict:
        """
        Parse MSA alignment file.
        :param file_path: path to alignment file
        :return: dictionary with ids as keys and sequences as values
        """

        def add_seq(aln: dict, sequence_id: str, seq_list: list):
            """
            Add a complete sequence and check for non-allowed chars
            :param aln: alignment dictionary to build
            :param sequence_id: sequence id to add
            :param seq_list: sequences to add
            :return: alignment with added sequences
            """
            final_seq = ''.join(seq_list).upper()
            # Check for non-allowed characters
            invalid_chars = set(final_seq) - set(config.POSSIBLE_CHARS)
            if invalid_chars:
                raise ValueError(
                    f"{sequence_id} contains invalid characters: {', '.join(invalid_chars)}. Allowed chars are: {config.POSSIBLE_CHARS}"
                )
            aln[sequence_id] = final_seq

            return aln

        alignment, seq_lines = {}, []
        seq_id = None

        with open(file_path, 'r') as file:
            for i, line in enumerate(file):
                line = line.strip()
                # initial check for fasta format
                if i == 0 and not line.startswith(">"):
                    raise ValueError('Alignment has to be in fasta format starting with >SeqID.')
                if line.startswith(">"):
                    if seq_id:
                        alignment = add_seq(alignment, seq_id, seq_lines)
                    # initialize a new sequence
                    seq_id, seq_lines = line[1:], []
                else:
                    seq_lines.append(line)
            # handle last sequence
            if seq_id:
                alignment = add_seq(alignment, seq_id, seq_lines)
        # final sanity checks
        if alignment:
            # alignment contains only one sequence:
            if len(alignment) < 2:
                raise ValueError("Alignment must contain more than one sequence.")
            # alignment sequences are not same length
            first_seq_len = len(next(iter(alignment.values())))
            for sequence_id, sequence in alignment.items():
                if len(sequence) != first_seq_len:
                    raise ValueError(
                        f"All alignment sequences must have the same length. Sequence '{sequence_id}' has length {len(sequence)}, expected {first_seq_len}."
                    )
            # all checks passed
            return alignment
        else:
            raise ValueError(f"Alignment file {file_path} does not contain any sequences in fasta format.")

    @staticmethod
    def _validate_ref(reference: str | None, alignment: dict) -> str | None | ValueError:
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
    def _validate_zoom(zoom: tuple | int, original_aln: dict) -> ValueError | tuple | None:
        """
        Validates if the user defined zoom range is within the start, end of the initial
        alignment.\n
        :param zoom: zoom range or zoom start
        :param original_aln: non-zoomed alignment dict
        :return: validated zoom range
        """
        if zoom is not None:
            aln_length = len(original_aln[list(original_aln.keys())[0]])
            # check if only over value is provided -> stop is alignment length
            if isinstance(zoom, int):
                if 0 <= zoom < aln_length:
                    return zoom, aln_length - 1
                else:
                    raise ValueError('Zoom start must be within the alignment length range.')
            # check if more than 2 values are provided
            if len(zoom) != 2:
                raise ValueError('Zoom position have to be (zoom_start, zoom_end)')
            # validate zoom start/stop
            for position in zoom:
                if type(position) != int:
                    raise ValueError('Zoom positions have to be integers.')
                if position not in range(0, aln_length):
                    raise ValueError('Zoom position out of range')

        return zoom

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
            if 'U' in alignment[record]:
                return 'RNA'
            counter += sum(map(alignment[record].count, ['A', 'C', 'G', 'T', 'N', '-']))
        # determine which is the most likely type
        if counter / len(alignment) >= 0.7 * len(alignment[list(alignment.keys())[0]]):
            return 'DNA'
        else:
            return 'AA'

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
    def zoom(self, zoom_pos: tuple | int):
        """
        Validate if the user defined zoom range.
        """
        self._zoom = self._validate_zoom(zoom_pos, self._alignment)

    # Property without setters
    @property
    def aln_type(self) -> str:
        """
        define the aln type:
        RNA, DNA or AA
        """
        return self._aln_type

    # On the fly properties without setters
    @property
    def length(self) -> int:
        return len(next(iter(self.alignment.values())))

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
    def get_reference_coords(self) -> tuple[int, int]:
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
                if self.alignment[self.reference_id][start] not in ['-', 'N']:
                    break
            # 3' --> 5'
            for end in range(self.length - 1, 0, -1):
                if self.alignment[self.reference_id][end] not in ['-', 'N']:
                    break

            return start, end

    def get_consensus(self, threshold: float = None, use_ambig_nt: bool = False) -> str:
        """
        Creates a non-gapped consensus sequence.

        threshold: Threshold for consensus sequence. If use_ambig_nt = True the ambig. char that encodes
        the nucleotides that reach a cumulative frequency >= threshold is used. Otherwise 'N' (for nt alignments)
        or 'X' (for as alignments) is used if none of the characters reach a cumulative frequency >= threshold.

        use_ambig_nt:test Use ambiguous character nt if none of the possible nt at a alignment position
        has a frequency above the defined threshold.

        :return: consensus sequence
        """

        # helper functions
        def determine_counts(alignment_dict: dict, position: int) -> dict:
            """
            count the number of each char at
            an idx of the alignment. return sorted dic.
            handles ambiguous nucleotides in sequences.
            also handles gaps.
            """
            nucleotide_list = []

            # get all nucleotides
            for sequence in alignment_dict.items():
                nucleotide_list.append(sequence[1][position])
            # count occurences of nucleotides
            counter = dict(collections.Counter(nucleotide_list))
            # get permutations of an ambiguous nucleotide
            to_delete = []
            temp_dict = {}
            for nucleotide in counter:
                if nucleotide in config.AMBIG_CHARS[self.aln_type]:
                    to_delete.append(nucleotide)
                    permutations = config.AMBIG_CHARS[self.aln_type][nucleotide]
                    adjusted_freq = 1 / len(permutations)
                    for permutation in permutations:
                        if permutation in temp_dict:
                            temp_dict[permutation] += adjusted_freq
                        else:
                            temp_dict[permutation] = adjusted_freq

            # drop ambiguous entries and add adjusted freqs to
            if to_delete:
                for i in to_delete:
                    counter.pop(i)
                for nucleotide in temp_dict:
                    if nucleotide in counter:
                        counter[nucleotide] += temp_dict[nucleotide]
                    else:
                        counter[nucleotide] = temp_dict[nucleotide]

            return dict(sorted(counter.items(), key=lambda x: x[1], reverse=True))

        def get_consensus_char(counts: dict, cutoff: float) -> list:
            """
            get a list of nucleotides for the consensus seq
            """
            n = 0

            consensus_chars = []
            for char in counts:
                n += counts[char]
                consensus_chars.append(char)
                if n >= cutoff:
                    break

            return consensus_chars

        def get_ambiguous_char(nucleotides: list) -> str:
            """
            get ambiguous char from a list of nucleotides
            """
            for ambiguous, permutations in config.AMBIG_CHARS[self.aln_type].items():
                if set(permutations) == set(nucleotides):
                    return ambiguous

        # check if params have been set correctly
        if threshold is not None:
            if threshold < 0 or threshold > 1:
                raise ValueError('Threshold must be between 0 and 1.')
        if self.aln_type == 'AA' and use_ambig_nt:
            raise ValueError('Ambiguous characters can not be calculated for amino acid alignments.')
        if threshold is None and use_ambig_nt:
            raise ValueError('To calculate ambiguous nucleotides, set a threshold > 0.')

        alignment = self.alignment
        consensus = str()

        if threshold is not None:
            consensus_cutoff = len(alignment) * threshold
        else:
            consensus_cutoff = 0

        # built consensus sequences
        for idx in range(self.length):
            char_counts = determine_counts(alignment, idx)
            consensus_chars = get_consensus_char(
                char_counts,
                consensus_cutoff
            )
            if threshold != 0:
                if len(consensus_chars) > 1:
                    if use_ambig_nt:
                        char = get_ambiguous_char(consensus_chars)
                    else:
                        if self.aln_type == 'AA':
                            char = 'X'
                        else:
                            char = 'N'
                    consensus = consensus + char
                else:
                    consensus = consensus + consensus_chars[0]
            else:
                consensus = consensus + consensus_chars[0]

        return consensus

    def get_conserved_orfs(self, min_length: int = 100, identity_cutoff: float | None = None) -> dict:
        """
        conserved ORF definition:
            - conserved starts and stops
            - start, stop must be on the same frame
            - stop - start must be at least min_length
            - all ungapped seqs[start:stop] must have at least min_length
            - no ungapped seq can have a Stop in between Start Stop

        Conservation is measured by number of positions with identical characters divided by
        orf slice of the alignment.

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

        # helper functions
        def determine_conserved_start_stops(alignment: dict, alignment_length: int) -> tuple:
            """
            Determine all start and stop codons within an alignment.
            :param alignment: alignment
            :param alignment_length: length of alignment
            :return: start and stop codon positions
            """
            starts = config.START_CODONS[self.aln_type]
            stops = config.STOP_CODONS[self.aln_type]

            list_of_starts, list_of_stops = [], []
            ref = alignment[list(alignment.keys())[0]]
            for nt_position in range(alignment_length):
                if ref[nt_position:nt_position + 3] in starts:
                    conserved_start = True
                    for sequence in alignment:
                        if not alignment[sequence][nt_position:].replace('-', '')[0:3] in starts:
                            conserved_start = False
                            break
                    if conserved_start:
                        list_of_starts.append(nt_position)

                if ref[nt_position:nt_position + 3] in stops:
                    conserved_stop = True
                    for sequence in alignment:
                        if not alignment[sequence][nt_position:].replace('-', '')[0:3] in stops:
                            conserved_stop = False
                            break
                    if conserved_stop:
                        list_of_stops.append(nt_position)

            return list_of_starts, list_of_stops

        def get_ungapped_sliced_seqs(alignment: dict, start_pos: int, stop_pos: int) -> list:
            """
            get ungapped sequences starting and stop codons and eliminate gaps
            :param alignment: alignment
            :param start_pos: start codon
            :param stop_pos: stop codon
            :return: sliced sequences
            """
            ungapped_seqs = []
            for seq_id in alignment:
                ungapped_seqs.append(alignment[seq_id][start_pos:stop_pos + 3].replace('-', ''))

            return ungapped_seqs

        def additional_stops(ungapped_seqs: list) -> bool:
            """
            Checks for the presence of a stop codon
            :param ungapped_seqs: list of ungapped sequences
            :return: Additional stop codons (True/False)
            """
            stops = config.STOP_CODONS[self.aln_type]

            for sliced_seq in ungapped_seqs:
                for position in range(0, len(sliced_seq) - 3, 3):
                    if sliced_seq[position:position + 3] in stops:
                        return True
            return False

        def calculate_identity(identity_matrix: ndarray, aln_slice:list) -> float:
            sliced_array = identity_matrix[:,aln_slice[0]:aln_slice[1]]
            return np.sum(np.all(sliced_array == 1, axis=0))/(aln_slice[1] - aln_slice[0]) * 100

        # checks for arguments

        if self.aln_type == 'AA':
            raise TypeError('ORF search only for RNA/DNA alignments')

        if identity_cutoff is not None:
            if identity_cutoff > 100 or identity_cutoff < 0:
                raise ValueError('conservation cutoff must be between 0 and 100')

        if min_length <= 0 or min_length > self.length:
            raise ValueError(f'min_length must be between 0 and {self.length}')

        # ini
        identities = self.calc_identity_alignment()
        alignments = [self.alignment, self.calc_reverse_complement_alignment()]
        aln_len = self.length

        orf_counter = 0
        orf_dict = {}

        for aln, direction in zip(alignments, ['+', '-']):
            # check for starts and stops in the first seq and then check if these are present in all seqs
            conserved_starts, conserved_stops = determine_conserved_start_stops(aln, aln_len)
            # check each frame
            for frame in (0, 1, 2):
                potential_starts = [x for x in conserved_starts if x % 3 == frame]
                potential_stops = [x for x in conserved_stops if x % 3 == frame]
                last_stop = -1
                for start in potential_starts:
                    # go to the next stop that is sufficiently far away in the alignment
                    next_stops = [x for x in potential_stops if x + 3 >= start + min_length]
                    if not next_stops:
                        continue
                    next_stop = next_stops[0]
                    ungapped_sliced_seqs = get_ungapped_sliced_seqs(aln, start, next_stop)
                    # re-check the lengths of all ungapped seqs
                    ungapped_seq_lengths = [len(x) >= min_length for x in ungapped_sliced_seqs]
                    if not all(ungapped_seq_lengths):
                        continue
                    # if no stop codon between start and stop --> write to dictionary
                    if not additional_stops(ungapped_sliced_seqs):
                        if direction == '+':
                            positions = [start, next_stop + 3]
                        else:
                            positions = [aln_len - next_stop - 3, aln_len - start]
                        if last_stop != next_stop:
                            last_stop = next_stop
                            conservation = calculate_identity(identities, positions)
                            if identity_cutoff is not None and conservation < identity_cutoff:
                                continue
                            orf_dict[f'ORF_{orf_counter}'] = {'location': [positions],
                                                              'frame': frame,
                                                              'strand': direction,
                                                              'conservation': conservation,
                                                              'internal': []
                                                              }
                            orf_counter += 1
                        else:
                            orf_dict[f'ORF_{orf_counter - 1}']['internal'].append(positions)

        return orf_dict

    def get_non_overlapping_conserved_orfs(self, min_length: int = 100, identity_cutoff:float = None) -> dict:
        """

        First calculates all ORFs and then searches from 5'
        all non-overlapping orfs in the fw strand and from the
        3' all non-overlapping orfs in th rw strand.

        No overlap algorithm:\n
        frame 1: -[M------*]--- ----[M--*]---------[M-----\n
        frame 2: -------[M------*]---------[M---*]--------\n
        frame 3: [M---*]-----[M----------*]----------[M---\n
        \n
        results: [M---*][M------*]--[M--*]-[M---*]-[M-----\n
        frame:    3      2           1      2       1

        :return: dictionary with non-overlapping orfs
        """
        orf_dict = self.get_conserved_orfs(min_length, identity_cutoff)

        fw_orfs, rw_orfs = [], []

        for orf in orf_dict:
            if orf_dict[orf]['strand'] == '+':
                fw_orfs.append((orf, orf_dict[orf]['location'][0]))
            else:
                rw_orfs.append((orf, orf_dict[orf]['location'][0]))

        fw_orfs.sort(key=lambda x: x[1][0])  # sort by start pos
        rw_orfs.sort(key=lambda x: x[1][1], reverse=True)  # sort by stop pos
        non_overlapping_orfs = []
        for orf_list, strand in zip([fw_orfs, rw_orfs], ['+', '-']):
            previous_stop = -1 if strand == '+' else self.length + 1
            for orf in orf_list:
                if strand == '+' and orf[1][0] > previous_stop:
                    non_overlapping_orfs.append(orf[0])
                    previous_stop = orf[1][1]
                elif strand == '-' and orf[1][1] < previous_stop:
                    non_overlapping_orfs.append(orf[0])
                    previous_stop = orf[1][0]

        non_overlap_dict = {}
        for orf in orf_dict:
            if orf in non_overlapping_orfs:
                non_overlap_dict[orf] = orf_dict[orf]

        return non_overlap_dict

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

        # helper functions
        def shannons_entropy(character_list: list, states: int, aln_type: str) -> float:
            """
            Calculate the shannon's entropy of a sequence and
            normalized between 0 and 1.
            :param character_list: characters at an alignment position
            :param states: number of potential characters that can be present
            :param aln_type: type of the alignment
            :returns: entropy
            """
            ent, n_chars = 0, len(character_list)
            # only one char is in the list
            if n_chars <= 1:
                return ent
            # calculate the number of unique chars and their counts
            chars, char_counts = np.unique(character_list, return_counts=True)
            char_counts = char_counts.astype(float)
            # ignore gaps for entropy calc
            char_counts, chars = char_counts[chars != "-"], chars[chars != "-"]
            # correctly handle ambiguous chars
            index_to_drop = []
            for index, char in enumerate(chars):
                if char in config.AMBIG_CHARS[aln_type]:
                    index_to_drop.append(index)
                    amb_chars, amb_counts = np.unique(config.AMBIG_CHARS[aln_type][char], return_counts=True)
                    amb_counts = amb_counts / len(config.AMBIG_CHARS[aln_type][char])
                    # add the proportionate numbers to initial array
                    for amb_char, amb_count in zip(amb_chars, amb_counts):
                        if amb_char in chars:
                            char_counts[chars == amb_char] += amb_count
                        else:
                            chars, char_counts = np.append(chars, amb_char), np.append(char_counts, amb_count)
            # drop the ambiguous characters from array
            char_counts, chars = np.delete(char_counts, index_to_drop), np.delete(chars, index_to_drop)
            # calc the entropy
            probs = char_counts / n_chars
            if np.count_nonzero(probs) <= 1:
                return ent
            for prob in probs:
                ent -= prob * math.log(prob, states)

            return ent

        aln = self.alignment
        entropys = []

        if self.aln_type == 'AA':
            states = 20
        else:
            states = 4
        # iterate over alignment positions and the sequences
        for nuc_pos in range(self.length):
            pos = []
            for record in aln:
                pos.append(aln[record][nuc_pos])
            entropys.append(shannons_entropy(pos, states, self.aln_type))

        return entropys

    def calc_gc(self) -> list | TypeError:
        """
        Determine the GC content for every position in an nt alignment.
        :return: GC content for every position.
        """
        if self.aln_type == 'AA':
            raise TypeError("GC computation is not possible for aminoacid alignment")

        gc, aln, amb_nucs = [], self.alignment, config.AMBIG_CHARS[self.aln_type]

        for position in range(self.length):
            nucleotides = str()
            for record in aln:
                nucleotides = nucleotides + aln[record][position]
            # ini dict with chars that occur and which ones to
            # count in which freq
            to_count = {
                'G': 1,
                'C': 1,
            }
            # handle ambig. nuc
            for char in amb_nucs:
                if char in nucleotides:
                    to_count[char] = (amb_nucs[char].count('C') + amb_nucs[char].count('G')) / len(amb_nucs[char])

            gc.append(
                sum([nucleotides.count(x) * to_count[x] for x in to_count]) / len(nucleotides)
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
            for record in aln.keys():
                pos = pos + aln[record][nuc_pos]
            coverage.append(1 - pos.count('-') / len(pos))

        return coverage

    def calc_reverse_complement_alignment(self) -> dict | TypeError:
        """
        Reverse complement the alignment.
        :return: Alignment
        """
        if self.aln_type == 'AA':
            raise TypeError('Reverse complement only for RNA or DNA.')

        aln = self.alignment
        reverse_complement_dict = {}

        for seq_id in aln:
            reverse_complement_dict[seq_id] = ''.join(config.COMPLEMENT[base] for base in reversed(aln[seq_id]))

        return reverse_complement_dict

    def calc_identity_alignment(self, encode_mismatches:bool=True, encode_mask:bool=False, encode_gaps:bool=True, encode_ambiguities:bool=False) -> np.ndarray:
        """
        Converts alignment to identity array (identical=1) compared to majority consensus or reference:\n

        :param encode_mismatches: encode mismatch as 0
        :param encode_mask: encode mask with value=2 --> also in the reference
        :param encode_gaps: encode gaps with np.nan --> also in the reference
        :param encode_ambiguities: encode ambiguities with value=3
        :return: identity alignment
        """

        aln = self.alignment
        ref = aln[self.reference_id] if self.reference_id is not None else self.get_consensus()

        # convert alignment to array
        sequences = np.array([list(aln[seq_id]) for seq_id in list(aln.keys())])
        reference = np.array(list(ref))
        # ini matrix
        identity_matrix = np.full(sequences.shape, 0, dtype=float)

        is_identical = sequences == reference

        if encode_gaps:
            is_gap = sequences == '-'
        else:
            is_gap = np.full(sequences.shape, False)

        if encode_mask:
            if self.aln_type == 'AA':
                is_n_or_x = np.isin(sequences, ['X'])
            else:
                is_n_or_x = np.isin(sequences, ['N'])
        else:
            is_n_or_x = np.full(sequences.shape, False)

        if encode_ambiguities:
            is_ambig = np.isin(sequences, [key for key in config.AMBIG_CHARS[self.aln_type] if key not in ['N', 'X', '-']])
        else:
            is_ambig = np.full(sequences.shape, False)

        if encode_mismatches:
            is_mismatch = ~is_gap & ~is_identical & ~is_n_or_x & ~is_ambig
        else:
            is_mismatch = np.full(sequences.shape, False)
        # assign values based on conditions
        identity_matrix[is_mismatch] = -1  # mismatch
        identity_matrix[is_gap] = np.nan  # gap
        identity_matrix[is_n_or_x] = -2  # 'N' or 'X'
        identity_matrix[is_ambig] = -3  # ambiguities

        return identity_matrix

    def calc_similarity_alignment(self, matrix_type:str|None=None) -> np.ndarray:
        """
        Calculate the similarity score between the alignment and the reference seq. Gaps are encoded as np.nan
        :param: matrix_type: type of similarity score (if not set - AA: BLOSSUM65, RNA/DNA: BLASTN)
        :return: identity array
        """

        aln = self.alignment
        ref = aln[self.reference_id] if self.reference_id is not None else self.get_consensus()
        if matrix_type is None:
            if self.aln_type == 'AA':
                matrix_type = 'BLOSUM65'
            else:
                matrix_type = 'TRANS'
        # load substitution matrix as dictionary
        try:
            subs_matrix = config.SUBS_MATRICES[self.aln_type][matrix_type]
        except KeyError:
            raise ValueError(
                f'The specified matrix does not exist for alignment type.\nAvailable matrices for {self.aln_type} are:\n{list(config.SUBS_MATRICES[self.aln_type].keys())}'
            )

        # set dtype and convert alignment to a NumPy array for vectorized processing
        dtype = np.dtype(float, metadata={'matrix': matrix_type})
        sequences = np.array([list(aln[seq_id]) for seq_id in list(aln.keys())])
        reference = np.array(list(ref))
        valid_chars = list(subs_matrix.keys())
        identity_array = np.full(sequences.shape, np.nan, dtype=dtype)

        for i, seq in enumerate(sequences):
            for j, char in enumerate(seq):
                if char in valid_chars and reference[j] in valid_chars:
                    identity_array[i, j] = subs_matrix[char][reference[j]]


        return identity_array

    def calc_percent_recovery(self) -> dict:
        """
        Recovery per sequence either compared to the majority consensus seq
        or the reference seq.\n
        Defined as:\n

        (1 - sum(N/X/- characters in ungapped ref regions))*100

        This is highly similar to how nextclade calculates recovery over reference.

        :return: dict
        """

        aln = self.alignment

        if self.reference_id is not None:
            ref = aln[self.reference_id]
        else:
            ref = self.get_consensus()  # majority consensus

        if not any(char != '-' for char in ref):
            raise ValueError("Reference sequence is entirely gapped, cannot calculate recovery.")


        # count 'N', 'X' and '-' chars in non-gapped regions
        recovery_over_ref = dict()

        # Get positions of non-gap characters in the reference
        non_gap_positions = [i for i, char in enumerate(ref) if char != '-']
        cumulative_length = len(non_gap_positions)

        # Calculate recovery
        for seq_id in aln:
            if seq_id == self.reference_id:
                continue
            seq = aln[seq_id]
            count_invalid = sum(
                seq[pos] == '-' or
                (seq[pos] == 'X' if self.aln_type == "AA" else seq[pos] == 'N')
                for pos in non_gap_positions
            )
            recovery_over_ref[seq_id] = (1 - count_invalid / cumulative_length) * 100

        return recovery_over_ref

    def calc_character_frequencies(self) -> dict:
        """
        Calculate the percentage characters in the alignment:
        The frequencies are counted by seq and in total. The
        percentage of non-gap characters in the alignment is
        relative to the total number of non-gap characters.
        The gap percentage is relative to the sequence length.

        The output is a nested dictionary.

        :return: Character frequencies
        """

        aln, aln_length = self.alignment, self.length

        freqs = {'total': {'-': {'counts': 0, '% of alignment': float()}}}

        for seq_id in aln:
            freqs[seq_id], all_chars = {'-': {'counts': 0, '% of alignment': float()}}, 0
            unique_chars = set(aln[seq_id])
            for char in unique_chars:
                if char == '-':
                    continue
                # add characters to dictionaries
                if char not in freqs[seq_id]:
                    freqs[seq_id][char] = {'counts': 0, '% of non-gapped': 0}
                if char not in freqs['total']:
                    freqs['total'][char] = {'counts': 0, '% of non-gapped': 0}
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

    def calc_pairwise_identity_matrix(self, distance_type:str='ghd') -> ndarray:
        """
        Calculate pairwise identities for an alignment. As there are different definitions of sequence identity, there are different options implemented:

        ghd (global hamming distance): At each alignment position, check if characters match:
        \ndistance = matches / alignment_length * 100

        lhd (local hamming distance): Restrict the alignment to the region in both sequences that do not start and end with gaps:
        \ndistance = matches / min(5'3' ungapped seq1, 5'3' ungapped seq2) * 100

        ged (gap excluded distance): All gaps are excluded from the alignment
        \ndistance = matches / (matches + mismatches) * 100

        gcd (gap compressed distance): All consecutive gaps are compressed to one mismatch.
        \ndistance = matches / gap_compressed_alignment_length * 100

        :return: array with pairwise distances.
        """

        def hamming_distance(seq1: str, seq2: str) -> int:
            return sum(c1 == c2 for c1, c2 in zip(seq1, seq2))

        def ghd(seq1: str, seq2: str) -> float:
            return hamming_distance(seq1, seq2) / self.length * 100

        def lhd(seq1: str, seq2: str) -> float:
            # remove 5' trailing gaps
            i, j = 0, self.length - 1
            while i < self.length and (seq1[i] == '-' or seq2[i] == '-'):
                i += 1
            while j >= 0 and (seq1[j] == '-' or seq2[j] == '-'):
                j -= 1
            if i > j:
                return 0.0
            # slice seq
            seq1_, seq2_ = seq1[i:j + 1], seq2[i:j + 1]

            return hamming_distance(seq1_, seq2_) / min([len(seq1_), len(seq2_)]) * 100

        def ged(seq1: str, seq2: str) -> float:

            matches, mismatches = 0, 0

            for c1, c2 in zip(seq1, seq2):
                if c1 != '-' and c2 != '-':
                    if c1 == c2:
                        matches += 1
                    else:
                        mismatches += 1
            return matches / (matches + mismatches) * 100 if (matches + mismatches) > 0 else 0

        def gcd(seq1: str, seq2: str) -> float:
            matches = 0
            mismatches = 0
            in_gap = False

            for char1, char2 in zip(seq1, seq2):
                if char1 == '-' and char2 == '-':  # Shared gap: do nothing
                    continue
                elif char1 == '-' or char2 == '-':  # Gap in only one sequence
                    if not in_gap:  # Start of a new gap stretch
                        mismatches += 1
                        in_gap = True
                else:  # No gaps
                    in_gap = False
                    if char1 == char2:  # Matching characters
                        matches += 1
                    else:  # Mismatched characters
                        mismatches += 1

            return matches / (matches + mismatches) * 100 if (matches + mismatches) > 0 else 0


        # Map distance type to corresponding function
        distance_functions: Dict[str, Callable[[str, str], float]] = {
            'ghd': ghd,
            'lhd': lhd,
            'ged': ged,
            'gcd': gcd
        }

        if distance_type not in distance_functions:
            raise ValueError(f"Invalid distance type '{distance_type}'. Choose from {list(distance_functions.keys())}.")

        # Compute pairwise distances
        aln = self.alignment
        distance_func = distance_functions[distance_type]
        distance_matrix = np.zeros((len(aln), len(aln)))

        sequences = list(aln.values())
        for i, seq1 in enumerate(sequences):
            for j, seq2 in enumerate(sequences):
                if i <= j:  # Compute only once for symmetric matrix
                    distance_matrix[i, j] = distance_func(seq1, seq2)
                    distance_matrix[j, i] = distance_matrix[i, j]

        return distance_matrix

    def get_snps(self, include_ambig:bool=False) -> dict:
        """
        Calculate snps similar to snp-sites (output is comparable):
        https://github.com/sanger-pathogens/snp-sites
        Importantly, SNPs are only considered if at least one of the snps is not an ambiguous character.
        The SNPs are compared to a majority consensus sequence or to a reference if it has been set.

        :param include_ambig: Include ambiguous snps (default: False)
        :return: dictionary containing snp positions and their variants including their frequency.
        """
        aln = self.alignment
        ref = aln[self.reference_id] if self.reference_id is not None else self.get_consensus()
        aln = {x: aln[x] for x in aln.keys() if x != self.reference_id}
        seq_ids = list(aln.keys())
        snp_dict = {'#CHROM': self.reference_id if self.reference_id is not None else 'consensus', 'POS': {}}

        for pos in range(self.length):
            reference_char = ref[pos]
            if not include_ambig:
                if reference_char in config.AMBIG_CHARS[self.aln_type] and reference_char != '-':
                    continue
            alt_chars, snps = [], []
            for i, seq_id in enumerate(aln.keys()):
                alt_chars.append(aln[seq_id][pos])
                if reference_char != aln[seq_id][pos]:
                    snps.append(i)
            if not snps:
                continue
            if include_ambig:
                if all(alt_chars[x] in config.AMBIG_CHARS[self.aln_type] for x in snps):
                    continue
            else:
                snps = [x for x in snps if alt_chars[x] not in config.AMBIG_CHARS[self.aln_type]]
                if not snps:
                    continue
            if pos not in snp_dict:
                snp_dict['POS'][pos] = {'ref': reference_char, 'ALT': {}}
            for snp in snps:
                if alt_chars[snp] not in snp_dict['POS'][pos]['ALT']:
                    snp_dict['POS'][pos]['ALT'][alt_chars[snp]] = {
                        'AF': 1,
                        'SEQ_ID': [seq_ids[snp]]
                    }
                else:
                    snp_dict['POS'][pos]['ALT'][alt_chars[snp]]['AF'] += 1
                    snp_dict['POS'][pos]['ALT'][alt_chars[snp]]['SEQ_ID'].append(seq_ids[snp])
            # calculate AF
            if pos in snp_dict['POS']:
                for alt in snp_dict['POS'][pos]['ALT']:
                    snp_dict['POS'][pos]['ALT'][alt]['AF'] /= len(aln)

        return snp_dict
    # TODO: add functions to manipulate alignments


class Annotation:
    """
    An annotation class that allows to read in gff, gb or bed files and adjust its locations to that of the MSA.
    """

    def __init__(self, aln: MSA, annotation_path: str):
        """
        The annotation class. Lets you parse multiple standard formats
        which might be used for annotating an alignment. The main purpose
        is to parse the annotation file and adapt the locations of diverse
        features to the locations within the alignment, considering the
        respective alignment positions. Importantly, IDs of the alignment
        and the MSA have to partly match.

        :param aln: MSA class
        :param annotation_path: path to annotation file (gb, bed, gff).

        """

        self.ann_type, self._seq_id, self.locus, self.features  = self._parse_annotation(annotation_path, aln)  # read annotation
        self._gapped_seq = self._MSA_validation_and_seq_extraction(aln, self._seq_id)  # extract gapped sequence
        self._position_map = self._build_position_map()  # build a position map
        self._map_to_alignment()  # adapt feature locations

    @staticmethod
    def _MSA_validation_and_seq_extraction(aln: MSA, seq_id: str) -> str:
        """
        extract gapped sequence from MSA that corresponds to annotation
        :param aln: MSA class
        :param seq_id: sequence id to extract
        :return: gapped sequence
        """
        if not isinstance(aln, MSA):
            raise ValueError('alignment has to be an MSA class. use explore.MSA() to read in alignment')
        else:
            return aln.alignment[seq_id]

    @staticmethod
    def _parse_annotation(annotation_path: str, aln: MSA) -> tuple[str, str, str, Dict]:

        def detect_annotation_type(file_path: str) -> str:
            """
            Detect the type of annotation file (GenBank, GFF, or BED) based
            on the first relevant line (excluding empty and #)

            :param file_path: Path to the annotation file.
            :return: The detected file type ('gb', 'gff', or 'bed').

            :raises ValueError: If the file type cannot be determined.
            """

            with open(file_path, 'r') as file:
                for line in file:
                    # skip empty lines and comments
                    if not line.strip() or line.startswith('#'):
                        continue
                   # genbank
                    if line.startswith('LOCUS'):
                        return 'gb'
                    # gff
                    if len(line.split('\t')) == 9:
                        # Check for expected values
                        columns = line.split('\t')
                        if columns[6] in ['+', '-', '.'] and re.match(r'^\d+$', columns[3]) and re.match(r'^\d+$',columns[4]):
                            return 'gff'
                    # BED files are tab-delimited with at least 3 fields: chrom, start, end
                    fields = line.split('\t')
                    if len(fields) >= 3 and re.match(r'^\d+$', fields[1]) and re.match(r'^\d+$', fields[2]):
                        return 'bed'
                    # only read in the first line
                    break

            raise ValueError(
                "File type could not be determined. Ensure the file follows a recognized format (GenBank, GFF, or BED).")

        def parse_gb(file_path) -> dict:
            """
            parse a genebank file to dictionary - primarily retained are the informations
            for qualifiers as these will be used for plotting.

            :param file_path: path to genebank file
            :return: nested dictionary

            """

            def sanitize_gb_location(string: str) -> tuple[list, str]:
                """
                see: https://www.insdc.org/submitting-standards/feature-table/
                """
                strand = '+'
                locations = []
                # check the direction of the annotation
                if 'complement' in string:
                    strand = '-'
                # sanitize operators
                for operator in ['complement(', 'join(', 'order(']:
                    string = string.strip(operator)
                # sanitize possible chars for splitting start stop -
                # however in the future might not simply do this
                # as some useful information is retained
                for char in ['>', '<', ')']:
                    string = string.replace(char, '')
                # check if we have multiple location e.g. due to splicing
                if ',' in string:
                    raw_locations = string.split(',')
                else:
                    raw_locations = [string]
                # try to split start and stop
                for location in raw_locations:
                    for sep in ['..', '.', '^']:
                        if sep in location:
                            sanitized_locations = [int(x) for x in location.split(sep)]
                            sanitized_locations[0] = sanitized_locations[0] - 1  # enforce 0-based starts
                            locations.append(sanitized_locations)
                            break

                return locations, strand


            records = {}
            with open(file_path, "r") as file:
                record = None
                in_features = False
                counter_dict = {}
                for line in file:
                    line = line.rstrip()
                    parts = line.split()
                    # extract the locus id
                    if line.startswith('LOCUS'):
                        if record:
                            records[record['locus']] = record
                        record = {
                            'locus': parts[1],
                            'features': {}
                        }

                    elif line.startswith('FEATURES'):
                        in_features = True

                    # ignore the sequence info
                    elif line.startswith('ORIGIN'):
                        in_features = False

                    # now write useful feature information to dictionary
                    elif in_features:
                        if not line.strip():
                            continue
                        if line[5] != ' ':
                            feature_type, qualifier = parts[0], parts[1]
                            if feature_type not in record['features']:
                                record['features'][feature_type] = {}
                                counter_dict[feature_type] = 0
                            locations, strand = sanitize_gb_location(qualifier)
                            record['features'][feature_type][counter_dict[feature_type]] = {
                                'location': locations,
                                'strand': strand
                            }
                            counter_dict[feature_type] += 1
                        else:
                            try:
                                qualifier_type, qualifier = parts[0].split('=')
                            except ValueError:  # we are in the coding sequence
                                qualifier = qualifier + parts[0]

                            qualifier_type, qualifier = qualifier_type.lstrip('/'), qualifier.strip('"')
                            last_index = counter_dict[feature_type] - 1
                            record['features'][feature_type][last_index][qualifier_type] = qualifier

            records[record['locus']] = record

            return records

        def parse_gff(file_path) -> dict:
            """
            Parse a GFF3 (General Feature Format) file into a dictionary structure.

            :param file_path: path to genebank file
            :return: nested dictionary

            """
            records = {}
            with open(file_path, 'r') as file:
                previous_id, previous_feature = None, None
                for line in file:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
                    # ensure that region and source features are not named differently for gff and gb
                    if feature_type == 'region':
                        feature_type = 'source'
                    if seqid not in records:
                        records[seqid] = {'locus': seqid, 'features': {}}
                    if feature_type not in records[seqid]['features']:
                        records[seqid]['features'][feature_type] = {}

                    feature_id = len(records[seqid]['features'][feature_type])
                    feature = {
                        'strand': strand,
                    }

                    # Parse attributes into key-value pairs
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            feature[key.strip()] = value.strip()

                    # check if feature are the same --> possible splicing
                    if previous_id is not None and previous_feature == feature:
                        records[seqid]['features'][feature_type][previous_id]['location'].append([int(start)-1, int(end)])
                    else:
                        records[seqid]['features'][feature_type][feature_id] = feature
                        records[seqid]['features'][feature_type][feature_id]['location'] = [[int(start) - 1, int(end)]]
                    # set new previous id and features -> new dict as 'location' is pointed in current feature and this
                    # is the only key different if next feature has the same entries
                    previous_id, previous_feature = feature_id, {key:value for key, value in feature.items() if key != 'location'}

            return records

        def parse_bed(file_path) -> dict:
            """
            Parse a BED file into a dictionary structure.

            :param file_path: path to genebank file
            :return: nested dictionary

            """
            records = {}
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    chrom, start, end, *optional = parts

                    if chrom not in records:
                        records[chrom] = {'locus': chrom, 'features': {}}
                    feature_type = 'region'
                    if feature_type not in records[chrom]['features']:
                        records[chrom]['features'][feature_type] = {}

                    feature_id = len(records[chrom]['features'][feature_type])
                    feature = {
                        'location': [[int(start), int(end)]],  # BED uses 0-based start, convert to 1-based
                        'strand': '+',  # assume '+' if not present
                    }

                    # Handle optional columns (name, score, strand) --> ignore 7-12
                    if len(optional) >= 1:
                        feature['name'] = optional[0]
                    if len(optional) >= 2:
                        feature['score'] = optional[1]
                    if len(optional) >= 3:
                        feature['strand'] = optional[2]

                    records[chrom]['features'][feature_type][feature_id] = feature

            return records

        parse_functions: Dict[str, Callable[[str], dict]] = {
            'gb': parse_gb,
            'bed': parse_bed,
            'gff': parse_gff,
        }
        # determine the annotation content -> should be standard formatted
        try:
            annotation_type = detect_annotation_type(annotation_path)
        except ValueError as err:
            raise err

        # read in the annotation
        annotations = parse_functions[annotation_type](annotation_path)

        # sanity check whether one of the annotation ids and alignment ids match
        annotation_found = False
        for annotation in annotations.keys():
            for aln_id in aln.alignment.keys():
                aln_id_sanitized = aln_id.split(' ')[0]
                # check in both directions
                if aln_id_sanitized in annotation:
                    annotation_found = True
                    break
                if annotation in aln_id_sanitized:
                    annotation_found = True
                    break

        if not annotation_found:
            raise ValueError(f'the annotations of {annotation_path} do not match any ids in the MSA')

        # return only the annotation that has been found, the respective type and the seq_id to map to
        return annotation_type, aln_id, annotations[annotation]['locus'], annotations[annotation]['features']


    def _build_position_map(self) -> Dict[int, int]:
        """
        build a position map from a sequence.

        :return genomic position: gapped position
        """

        position_map = {}
        genomic_pos = 0
        for aln_index, char in enumerate(self._gapped_seq):
            if char != '-':
                position_map[genomic_pos] = aln_index
                genomic_pos += 1
        # ensure the last genomic position is included
        position_map[genomic_pos] = len(self._gapped_seq)

        return position_map


    def _map_to_alignment(self):
        """
        Adjust all feature locations to alignment positions
        """

        def map_location(position_map: Dict[int, int], locations: list) -> list:
            """
            Map genomic locations to alignment positions using a precomputed position map.

            :param position_map: Positions mapped from gapped to ungapped
            :param locations: List of genomic start and end positions.
            :return: List of adjusted alignment positions.
            """

            aligned_locs = []
            for start, end in locations:
                try:
                    aligned_start = position_map[start]
                    aligned_end = position_map[end]
                    aligned_locs.append([aligned_start, aligned_end])
                except KeyError:
                    raise ValueError(f"Positions {start}-{end} lie outside of the position map.")

            return aligned_locs

        for feature_type, features in self.features.items():
            for feature_id, feature_data in features.items():
                original_locations = feature_data['location']
                aligned_locations = map_location(self._position_map, original_locations)
                feature_data['location'] = aligned_locations

    # TODO: transfer of annotations to other sequences (return respective ungapped seq and adapted annotation)