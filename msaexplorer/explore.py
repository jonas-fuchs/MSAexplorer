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
from typing import Callable, Dict

# installed
import numpy as np
from numpy import ndarray
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqIO.InsdcIO import GenBankIterator

# msaexplorer
from msaexplorer import config
from msaexplorer._data_classes import PairwiseDistance, AlignmentStats, LengthStats, OpenReadingFrame, OrfCollection, SingleNucleotidePolymorphism, VariantCollection
from msaexplorer._helpers import _get_line_iterator, _create_distance_calculation_function_mapping, _read_alignment

#TODO: Move outputs to dataclasses
class MSA:
    """
    An alignment class that allows computation of several stats. Supported inputs are file paths to alignments in "fasta",
    "clustal", "phylip", "stockholm", "nexus" formats, raw alignment strings, or Bio.Align.MultipleSeqAlignment for
    compatibility with Biopython.
    """

    # dunder methods
    def __init__(self, alignment_string: str | MultipleSeqAlignment, reference_id: str = None, zoom_range: tuple | int = None):
        """
        Initialize an Alignment object.
        :param alignment_string: Path to alignment file or raw alignment string
        :param reference_id: reference id
        :param zoom_range: start and stop positions to zoom into the alignment
        """
        self._alignment = _read_alignment(alignment_string)
        self._reference_id = self._validate_ref(reference_id, self._alignment)
        self._zoom = self._validate_zoom(zoom_range, self._alignment)
        self._aln_type = self._determine_aln_type(self._alignment)

    def __len__(self) -> int:
        return len(self.sequence_ids)

    def __getitem__(self, key: str) -> str:
        return self.alignment[key]

    def __contains__(self, item: str) -> bool:
        return item in self.sequence_ids

    def __iter__(self):
        return iter(self.alignment)

    # Static methods
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
    def reference_id(self) -> str:
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

    @property
    def sequence_ids(self) -> list:
        return list(self.alignment.keys())

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

    def _create_position_stat_result(self, stat_name: str, values: list | ndarray) -> AlignmentStats:
        """
        Build a shared dataclass for position-wise statistics.
        """
        values_array = np.asarray(values, dtype=float)
        if self.zoom is None:
            positions = np.arange(self.length, dtype=int)
        else:
            positions = np.arange(self.zoom[0], self.zoom[0] + self.length, dtype=int)

        return AlignmentStats(
            stat_name=stat_name,
            positions=positions,
            values=values_array,
        )

    def _to_array(self) -> ndarray:
        """convert alignment to a numpy array"""
        return np.array([list(self.alignment[seq_id]) for seq_id in self.sequence_ids])

    def _get_reference_seq(self) -> str:
        """get the sequence of the reference sequence or majority consensus"""
        return self.alignment[self.reference_id] if self.reference_id is not None else self.get_consensus()

    def get_reference_coords(self) -> tuple[int, int]:
        """
        Determine the start and end coordinates of the reference sequence
        defined as the first/last nucleotide in the reference sequence
        (excluding N and gaps).

        :return: Start, End
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

        :param threshold: Threshold for consensus sequence. If use_ambig_nt = True the ambig. char that encodes
            the nucleotides that reach a cumulative frequency >= threshold is used. Otherwise 'N' (for nt alignments)
            or 'X' (for as alignments) is used if none of the characters reach a cumulative frequency >= threshold.
        :param use_ambig_nt: Use ambiguous character nt if none of the possible nt at a alignment position
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
                    break

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

    def get_conserved_orfs(self, min_length: int = 100, identity_cutoff: float | None = None) -> OrfCollection:
        """
        **conserved ORF definition:**
            - conserved starts and stops
            - start, stop must be on the same frame
            - stop - start must be at least min_length
            - all ungapped seqs[start:stop] must have at least min_length
            - no ungapped seq can have a Stop in between Start Stop

        Conservation is measured by the number of positions with identical characters divided by
        orf slice of the alignment.

        **Algorithm overview:**
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
            # define one sequence (first) as reference (it does not matter which one)
            ref = alignment[self.sequence_ids[0]]
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
            sliced_array = identity_matrix[:,aln_slice[0]:aln_slice[1]] + 1  # identical = 0, different = -1 --> add 1
            return np.sum(np.all(sliced_array == 1, axis=0))/(aln_slice[1] - aln_slice[0]) * 100

        # checks for arguments
        if self.aln_type == 'AA':
            raise TypeError('ORF search only for RNA/DNA alignments')

        if identity_cutoff is not None:
            if identity_cutoff > 100 or identity_cutoff < 0:
                raise ValueError('conservation cutoff must be between 0 and 100')

        if min_length <= 6 or min_length > self.length:
            raise ValueError(f'min_length must be between 6 and {self.length}')

        # ini
        identities = self.calc_identity_alignment()
        alignments = [self.alignment, self.calc_reverse_complement_alignment()]
        aln_len = self.length

        orf_counter = 0
        # use mutable dicts during construction and convert to dataclass at the end for immutability
        temp_orfs: list[dict] = []

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
                            positions = (start, next_stop + 3)
                        else:
                            positions = (aln_len - next_stop - 3, aln_len - start)
                        if last_stop != next_stop:
                            last_stop = next_stop
                            conservation = calculate_identity(identities, list(positions))
                            if identity_cutoff is not None and conservation < identity_cutoff:
                                continue
                            temp_orfs.append({
                                'orf_id': f'ORF_{orf_counter}',
                                'location': [positions],
                                'frame': frame,
                                'strand': direction,
                                'conservation': conservation,
                                'internal': [],
                            })
                            orf_counter += 1
                        else:
                            if temp_orfs:
                                temp_orfs[-1]['internal'].append(positions)

        # convert mutable intermediate dicts to frozen dataclasses
        orf_list = [
            OpenReadingFrame(
                orf_id=t['orf_id'],
                location=tuple(t['location']),
                frame=t['frame'],
                strand=t['strand'],
                conservation=t['conservation'],
                internal=tuple(t['internal']),
            )
            for t in temp_orfs
        ]
        return OrfCollection(orfs=tuple(orf_list))

    def get_non_overlapping_conserved_orfs(self, min_length: int = 100, identity_cutoff:float = None) -> OrfCollection:
        """
        First calculates all ORFs and then searches from 5'
        all non-overlapping orfs in the fw strand and from the
        3' all non-overlapping orfs in th rw strand.

        **No overlap algorithm:**
            **frame 1:** -[M------*]--- ----[M--*]---------[M-----

            **frame 2:** -------[M------*]---------[M---*]--------

            **frame 3:** [M---*]-----[M----------*]----------[M---

            **results:** [M---*][M------*]--[M--*]-[M---*]-[M-----

            frame:    3      2           1      2       1

        :return: OrfContainer with non-overlapping orfs
        """
        all_orfs = self.get_conserved_orfs(min_length, identity_cutoff)

        fw_orfs, rw_orfs = [], []

        for orf in all_orfs:
            orf_obj = all_orfs[orf]
            if orf_obj.strand == '+':
                fw_orfs.append((orf, orf_obj.location[0]))
            else:
                rw_orfs.append((orf, orf_obj.location[0]))

        fw_orfs.sort(key=lambda x: x[1][0])  # sort by start pos
        rw_orfs.sort(key=lambda x: x[1][1], reverse=True)  # sort by stop pos
        non_overlapping_ids = []
        for orf_list, strand in zip([fw_orfs, rw_orfs], ['+', '-']):
            previous_stop = -1 if strand == '+' else self.length + 1
            for orf in orf_list:
                if strand == '+' and orf[1][0] >= previous_stop:
                    non_overlapping_ids.append(orf[0])
                    previous_stop = orf[1][1]
                elif strand == '-' and orf[1][1] <= previous_stop:
                    non_overlapping_ids.append(orf[0])
                    previous_stop = orf[1][0]

        return OrfCollection(orfs=tuple(all_orfs[orf_id] for orf_id in all_orfs if orf_id in non_overlapping_ids))

    def calc_length_stats(self) -> LengthStats:
        """
        Determine the stats for the length of the ungapped seqs in the alignment.
        :return: dataclass with length stats
        """

        seq_lengths = [len(self.alignment[x].replace('-', '')) for x in self.alignment]

        return LengthStats(
            n_sequences=len(self.alignment),
            mean_length=float(np.mean(seq_lengths)),
            std_length=float(np.std(seq_lengths)),
            min_length=int(np.min(seq_lengths)),
            max_length=int(np.max(seq_lengths)),
        )

    def calc_entropy(self) -> AlignmentStats:
        """
        Calculate the normalized shannon's entropy for every position in an alignment:

        - 1: high entropy
        - 0: low entropy

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

        return self._create_position_stat_result('entropy', entropys)

    def calc_gc(self) -> AlignmentStats | TypeError:
        """
        Determine the GC content for every position in an nt alignment.
        :return: GC content for every position.
        :raises: TypeError for AA alignments
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

        return self._create_position_stat_result('gc', gc)

    def calc_coverage(self) -> AlignmentStats:
        """
        Determine the coverage of every position in an alignment.
        This is defined as:
            1 - cumulative length of '-' characters

        :return: Coverage at each alignment position.
        """
        coverage, aln = [], self.alignment

        for nuc_pos in range(self.length):
            pos = str()
            for record in self.sequence_ids:
                pos = pos + aln[record][nuc_pos]
            coverage.append(1 - pos.count('-') / len(pos))

        return self._create_position_stat_result('coverage', coverage)

    def calc_gap_frequency(self) -> AlignmentStats:
        """
        Determine the gap frequency for every position in an alignment. This is the inverted coverage.
        """
        coverage = self.calc_coverage()

        return self._create_position_stat_result('gap frequency', 1 - coverage.values)

    def calc_reverse_complement_alignment(self) -> dict | TypeError:
        """
        Reverse complement the alignment.
        :return: Alignment (rv)
        """
        if self.aln_type == 'AA':
            raise TypeError('Reverse complement only for RNA or DNA.')

        aln = self.alignment
        reverse_complement_dict = {}

        for seq_id in aln:
            reverse_complement_dict[seq_id] = ''.join(config.COMPLEMENT[base] for base in reversed(aln[seq_id]))

        return reverse_complement_dict

    def calc_numerical_alignment(self, encode_mask:bool=False, encode_ambiguities:bool=False) -> ndarray:
        """
        Transforms the alignment to numerical values. Ambiguities are encoded as -3, mask as -2 and the
        remaining chars with the idx + 1 of config.CHAR_COLORS[self.aln_type]['standard'].

        :param encode_ambiguities: encode ambiguities as -2
        :param encode_mask: encode mask with as -3
        :returns matrix
        """

        sequences = self._to_array()
        # ini matrix
        numerical_matrix = np.full(sequences.shape, np.nan, dtype=float)
        # first encode mask
        if encode_mask:
            if self.aln_type == 'AA':
                is_n_or_x = np.isin(sequences, ['X'])
            else:
                is_n_or_x = np.isin(sequences, ['N'])
            numerical_matrix[is_n_or_x] = -2
        # next encode ambig chars
        if encode_ambiguities:
            numerical_matrix[np.isin(sequences, [key for key in config.AMBIG_CHARS[self.aln_type] if key not in ['N', 'X', '-']])] = -3
        # next convert each char into their respective values
        for idx, char in enumerate(config.CHAR_COLORS[self.aln_type]['standard']):
            numerical_matrix[np.isin(sequences, [char])] = idx + 1

        return numerical_matrix

    def calc_identity_alignment(self, encode_mismatches:bool=True, encode_mask:bool=False, encode_gaps:bool=True, encode_ambiguities:bool=False, encode_each_mismatch_char:bool=False) -> ndarray:
        """
        Converts alignment to identity array (identical=0) compared to majority consensus or reference:\n

        :param encode_mismatches: encode mismatch as -1
        :param encode_mask: encode mask with value=-2 --> also in the reference
        :param encode_gaps: encode gaps with np.nan --> also in the reference
        :param encode_ambiguities: encode ambiguities with value=-3
        :param encode_each_mismatch_char: for each mismatch encode characters separately - these values represent the idx+1 values of config.CHAR_COLORS[self.aln_type]['standard']
        :return: identity alignment
        """

        ref = self._get_reference_seq()

        # convert alignment to array
        sequences = self._to_array()
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

        # encode every different character
        if encode_each_mismatch_char:
            for idx, char in enumerate(config.CHAR_COLORS[self.aln_type]['standard']):
                new_encoding = np.isin(sequences, [char]) & is_mismatch
                identity_matrix[new_encoding] = idx + 1
        # or encode different with a single value
        else:
            identity_matrix[is_mismatch] = -1  # mismatch

        identity_matrix[is_gap] = np.nan  # gap
        identity_matrix[is_n_or_x] = -2  # 'N' or 'X'
        identity_matrix[is_ambig] = -3  # ambiguities

        return identity_matrix

    def calc_similarity_alignment(self, matrix_type:str|None=None, normalize:bool=True) -> ndarray:
        """
        Calculate the similarity score between the alignment and the reference sequence, with normalization to highlight
        differences. The similarity scores are scaled to the range [0, 1] based on the substitution matrix values for the
        reference residue at each column. Gaps are encoded as np.nan.

        The calculation follows these steps:

        1. **Reference Sequence**: If a reference sequence is provided (via `self.reference_id`), it is used. Otherwise,
           a consensus sequence is generated to serve as the reference.
        2. **Substitution Matrix**: The similarity between residues is determined using a substitution matrix, such as
           BLOSUM65 for amino acids or BLASTN for nucleotides. The matrix is loaded based on the alignment type.
        3. **Per-Column Normalization (optional)**:

        For each column in the alignment:
            - The residue in the reference sequence is treated as the baseline for that column.
            - The substitution scores for the reference residue are extracted from the substitution matrix.
            - The scores are normalized to the range [0, 1] using the minimum and maximum possible scores for the reference residue.
            - This ensures that identical residues (or those with high similarity to the reference) have high scores,
            while more dissimilar residues have lower scores.
        4. **Output**:

           - The normalized similarity scores are stored in a NumPy array.
           - Gaps (if any) or residues not present in the substitution matrix are encoded as `np.nan`.

        :param: matrix_type: type of similarity score (if not set - AA: BLOSSUM65, RNA/DNA: BLASTN)
        :param: normalize: whether to normalize the similarity scores to range [0, 1]
        :return: A 2D NumPy array where each entry corresponds to the normalized similarity score between the aligned residue
            and the reference residue for that column. Values range from 0 (low similarity) to 1 (high similarity).
            Gaps and invalid residues are encoded as `np.nan`.
        :raise: ValueError
            If the specified substitution matrix is not available for the given alignment type.
        """

        ref = self._get_reference_seq()

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
        sequences = self._to_array()
        reference = np.array(list(ref))
        valid_chars = list(subs_matrix.keys())
        similarity_array = np.full(sequences.shape, np.nan, dtype=dtype)

        for j, ref_char in enumerate(reference):
            if ref_char not in valid_chars + ['-']:
                continue
            # Get local min and max for the reference residue
            if normalize and ref_char != '-':
                local_scores = subs_matrix[ref_char].values()
                local_min, local_max = min(local_scores), max(local_scores)

            for i, char in enumerate(sequences[:, j]):
                if char not in valid_chars:
                    continue
                # classify the similarity as max if the reference has a gap
                similarity_score = subs_matrix[char][ref_char] if ref_char != '-' else 1
                similarity_array[i, j] = (similarity_score - local_min) / (local_max - local_min) if normalize and ref_char != '-' else similarity_score

        return similarity_array

    def calc_position_matrix(self, matrix_type:str='PWM') -> None | ndarray | ValueError:
        """
        Calculates a position matrix of the specified type for the given alignment. The function
        supports generating matrices of types Position Frequency Matrix (PFM), Position Probability
        Matrix (PPM), Position Weight Matrix (PWM), and cumulative Information Content (IC). It validates
        the provided matrix type and includes pseudo-count adjustments to ensure robust calculations.

        :param matrix_type: Type of position matrix to calculate. Accepted values are 'PFM', 'PPM',
            'PWM', and 'IC'. Defaults to 'PWM'.
        :type matrix_type: str
        :raises ValueError: If the provided `matrix_type` is not one of the accepted values.
        :return: A numpy array representing the calculated position matrix of the specified type.
        :rtype: np.ndarray
        """

        # ini
        if matrix_type not in ['PFM', 'PPM', 'IC', 'PWM']:
            raise ValueError('Matrix_type must be PFM, PPM, IC or PWM.')
        possible_chars = list(config.CHAR_COLORS[self.aln_type]['standard'].keys())
        sequences = self._to_array()

        # calc position frequency matrix
        pfm = np.array([np.sum(sequences == char, 0) for char in possible_chars])
        if matrix_type == 'PFM':
            return pfm

        # calc position probability matrix (probability)
        pseudo_count = 0.0001  # to avoid 0 values
        pfm = pfm + pseudo_count
        ppm_non_char_excluded = pfm/np.sum(pfm, axis=0)  # use this for pwm/ic calculation
        ppm = pfm/len(self.sequence_ids)  # calculate the frequency based on row number
        if matrix_type == 'PPM':
            return ppm

        # calc position weight matrix (log-likelihood)
        pwm = np.log2(ppm_non_char_excluded * len(possible_chars))
        if matrix_type == 'PWM':
            return pwm

        # calc information content per position (in bits) - can be used to scale a ppm for sequence logos
        ic = np.sum(ppm_non_char_excluded * pwm, axis=0)
        if matrix_type == 'IC':
            return ic

    def calc_percent_recovery(self) -> dict[str, float]:
        """
        Recovery per sequence either compared to the majority consensus seq
        or the reference seq.\n
        Defined as:\n

        `(1 - sum(N/X/- characters in ungapped ref regions))*100`

        This is highly similar to how nextclade calculates recovery over reference.

        :return: dict
        """

        aln = self.alignment
        ref = self._get_reference_seq()

        if not any(char != '-' for char in ref):
            raise ValueError("Reference sequence is entirely gapped, cannot calculate recovery.")


        # count 'N', 'X' and '-' chars in non-gapped regions
        recovery_over_ref = dict()

        # Get positions of non-gap characters in the reference
        non_gap_positions = [i for i, char in enumerate(ref) if char != '-']
        cumulative_length = len(non_gap_positions)

        # Calculate recovery
        for seq_id in self.sequence_ids:
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

    def calc_pairwise_identity_matrix(self, distance_type:str='ghd') -> PairwiseDistance:
        """
        Calculate pairwise identities for an alignment. Different options are implemented:

        **1) ghd (global hamming distance)**: At each alignment position, check if characters match:
        \ndistance = matches / alignment_length * 100

        **2) lhd (local hamming distance)**: Restrict the alignment to the region in both sequences that do not start and end with gaps:
        \ndistance = matches / min(5'3' ungapped seq1, 5'3' ungapped seq2) * 100

        **3) ged (gap excluded distance)**: All gaps are excluded from the alignment
        \ndistance = (1 - mismatches / total) * 100

        **4) gcd (gap compressed distance)**: All consecutive gaps are compressed to one mismatch.
        \ndistance = matches / gap_compressed_alignment_length * 100

        RNA/DNA only:

        **5) jc69 (Jukes-Cantor 1969)**: Gaps excluded. Applies the JC69 substitution model to correct
        the p-distance for multiple hits (assumes equal base frequencies and substitution rates).
        \ncorrected_identity = (1 - d_JC69) * 100,  where d = -(3/4) * ln(1 - (4/3) * p)

        **6) k2p (Kimura 2-Parameter / K80)**: Gaps excluded. Distinguishes transitions (Ts) and
        transversions (Tv). Returns (1 - d_K2P) * 100 as corrected percent identity.
        \nd = -(1/2) * ln(1 - 2P - Q) - (1/4) * ln(1 - 2Q),  P = Ts/total, Q = Tv/total

        :param distance_type: type of distance computation: ghd, lhd, ged, gcd and nucleotide only: jc69 and k2p
        :return: array with pairwise distances.
        """

        distance_functions = _create_distance_calculation_function_mapping()

        if distance_type not in distance_functions:
            raise ValueError(f"Invalid distance type '{distance_type}'. Choose from {list(distance_functions.keys())}.")

        aln = self.alignment

        if self.aln_type == 'AA' and distance_type in ['jc69', 'k2p']:
            raise ValueError(f"JC69 and K2P are not supported for {self.aln_type} alignment.")

        # Compute pairwise distances
        distance_func = distance_functions[distance_type]
        distance_matrix = np.zeros((len(aln), len(aln)))

        sequences = list(aln.values())
        n = len(sequences)
        for i in range(n):
            seq1 = sequences[i]
            for j in range(i, n):
                seq2 = sequences[j]
                dist = distance_func(seq1, seq2, self.length)
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist

        return PairwiseDistance(
            reference_id=None,
            sequence_ids=self.sequence_ids,
            distances=distance_matrix
        )

    def calc_pairwise_distance_to_reference(self, distance_type:str='ghd') -> PairwiseDistance:
        """
        Calculate pairwise identities between reference and all sequences in the alignment. Same computation as calc_pairwise_identity_matrix but compared to a single sequence. Different options are implemented:

        **1) ghd (global hamming distance)**: At each alignment position, check if characters match:
        \ndistance = matches / alignment_length * 100

        **2) lhd (local hamming distance)**: Restrict the alignment to the region in both sequences that do not start and end with gaps:
        \ndistance = matches / min(5'3' ungapped seq1, 5'3' ungapped seq2) * 100

        **3) ged (gap excluded distance)**: All gaps are excluded from the alignment
        \ndistance = (1 - mismatches / total) * 100

        **4) gcd (gap compressed distance)**: All consecutive gaps are compressed to one mismatch.
        \ndistance = matches / gap_compressed_alignment_length * 100

        RNA/DNA only:

        **5) jc69 (Jukes-Cantor 1969)**: Gaps excluded. Applies the JC69 substitution model to correct
        the p-distance for multiple hits (assumes equal base frequencies and substitution rates).
        \ncorrected_identity = (1 - d_JC69) * 100,  where d = -(3/4) * ln(1 - (4/3) * p)

        **6) k2p (Kimura 2-Parameter / K80)**: Gaps excluded. Distinguishes transitions (Ts) and
        transversions (Tv). Returns (1 - d_K2P) * 100 as corrected percent identity.
        \nd = -(1/2) * ln(1 - 2P - Q) - (1/4) * ln(1 - 2Q),  P = Ts/total, Q = Tv/total

        :param distance_type: type of distance computation: ghd, lhd, ged, gcd and nucleotide only: jc69 and k2p
        :return: array with pairwise distances.
        :return: dataclass with reference label, sequence ids and pairwise distances.
        """

        distance_functions = _create_distance_calculation_function_mapping()

        if distance_type not in distance_functions:
            raise ValueError(f"Invalid distance type '{distance_type}'. Choose from {list(distance_functions.keys())}.")

        if self.aln_type == 'AA' and distance_type in ['jc69', 'k2p']:
            raise ValueError(f"JC69 and K2P are not supported for {self.aln_type} alignment.")

        # Compute pairwise distances
        distance_func = distance_functions[distance_type]
        aln = self.alignment
        ref_id = self.reference_id

        ref_seq = aln[ref_id] if ref_id is not None else self.get_consensus()
        distances = []
        distance_names = []

        for seq_id in aln:
            if seq_id == ref_id:
                continue
            distance_names.append(seq_id)
            distances.append(distance_func(ref_seq, aln[seq_id], self.length))

        return PairwiseDistance(
            reference_id=ref_id if ref_id is not None else 'consensus',
            sequence_ids=distance_names,
            distances=np.array(distances)
        )

    def get_snps(self, include_ambig:bool=False) -> VariantCollection:
        """
        Calculate snps similar to snp-sites (output is comparable):
        https://github.com/sanger-pathogens/snp-sites
        Importantly, SNPs are only considered if at least one of the snps is not an ambiguous character.
        The SNPs are compared to a majority consensus sequence or to a reference if it has been set.

        :param include_ambig: Include ambiguous snps (default: False)
        :return: dataclass containing SNP positions and their variants including frequency.
        """
        aln = self.alignment
        ref = self._get_reference_seq()
        aln = {x: aln[x] for x in self.sequence_ids if x != self.reference_id}
        seq_ids = list(aln.keys())
        chrom = self.reference_id if self.reference_id is not None else 'consensus'
        snp_positions = {}
        aln_size = len(aln)

        for pos in range(self.length):
            reference_char = ref[pos]
            if not include_ambig:
                if reference_char in config.AMBIG_CHARS[self.aln_type] and reference_char != '-':
                    continue
            alt_chars, snps = [], []
            for i, seq_id in enumerate(seq_ids):
                alt_char = aln[seq_id][pos]
                alt_chars.append(alt_char)
                if reference_char != alt_char:
                    snps.append(i)
            if not snps:
                continue
            # Filter out ambiguous snps if not included
            if include_ambig:
                if all(alt_chars[x] in config.AMBIG_CHARS[self.aln_type] for x in snps):
                    continue
            else:
                snps = [x for x in snps if alt_chars[x] not in config.AMBIG_CHARS[self.aln_type]]
                if not snps:
                    continue
            # Build allele dict with counts
            alt_dict = {}
            for snp_idx in snps:
                alt_char = alt_chars[snp_idx]
                if alt_char not in alt_dict:
                    alt_dict[alt_char] = {'count': 0, 'seq_ids': []}
                alt_dict[alt_char]['count'] += 1
                alt_dict[alt_char]['seq_ids'].append(seq_ids[snp_idx])
            
            # Convert to final format: alt_char -> (frequency, seq_ids_tuple)
            alt = {
                alt_char: (data['count'] / aln_size, tuple(data['seq_ids']))
                for alt_char, data in alt_dict.items()
            }
            snp_positions[pos] = SingleNucleotidePolymorphism(ref=reference_char, alt=alt)

        return VariantCollection(chrom=chrom, positions=snp_positions)

    def calc_transition_transversion_score(self) -> AlignmentStats:
        """
        Based on the snp positions, calculates a transition/transversions score.
        A positive score means higher ratio of transitions and negative score means
        a higher ratio of transversions.
        :return: list
        """

        if self.aln_type == 'AA':
            raise TypeError('TS/TV scoring only for RNA/DNA alignments')

        # ini
        snps = self.get_snps()
        score = [0]*self.length

        for pos, snp in snps.positions.items():
            for alt, (af, _) in snp.alt.items():
                # check the type of substitution
                if snp.ref + alt in ['AG', 'GA', 'CT', 'TC', 'CU', 'UC']:
                    score[pos] += af
                else:
                    score[pos] -= af

        return self._create_position_stat_result('ts tv score', score)


class Annotation:
    """
    An annotation class that allows to read in gff, gb, or bed files and adjust its locations to that of the MSA.
    """

    def __init__(self, aln: MSA, annotation: str | GenBankIterator):
        """
        The annotation class. Lets you parse multiple standard formats
        which might be used for annotating an alignment. The main purpose
        is to parse the annotation file and adapt the locations of diverse
        features to the locations within the alignment, considering the
        respective alignment positions. Importantly, IDs of the alignment
        and the MSA have to partly match.

        :param aln: MSA class
        :param annotation: path to file (gb, bed, gff) or raw string or GenBankIterator from biopython

        """

        self.ann_type, self._seq_id, self.locus, self.features  = self._parse_annotation(annotation, aln)  # read annotation
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
            return aln._alignment[seq_id]

    @staticmethod
    def _parse_annotation(annotation: str | GenBankIterator, aln: MSA) -> tuple[str, str, str, Dict]:

        def detect_annotation_type(handle: str | GenBankIterator) -> str:
            """
            Detect the type of annotation file (GenBank, GFF, or BED) based
            on the first relevant line (excluding empty and #). Also recognizes
            Bio.SeqIO iterators as GenBank format.

            :param handle: Path to the annotation file or Bio.SeqIO iterator for genbank records read with biopython.
            :return: The detected file type ('gb', 'gff', or 'bed').

            :raises ValueError: If the file type cannot be determined.
            """
            # Check if input is a SeqIO iterator
            if isinstance(handle, GenBankIterator):
                return 'gb'

            with _get_line_iterator(handle) as h:
                for line in h:
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

        def parse_gb(file: str | GenBankIterator) -> dict:
            """
            Parse a GenBank file into the same dictionary structure used by the annotation pipeline.

            :param file: path to genbank file, raw string, or Bio.SeqIO iterator
            :return: nested dictionary
            """
            records = {}

            # Check if input is a GenBankIterator
            if isinstance(file, GenBankIterator):
                # Direct GenBankIterator input
                seq_records = list(file)
            else:
                # File path or string input
                with _get_line_iterator(file) as handle:
                    seq_records = list(SeqIO.parse(handle, "genbank"))

            for seq_record in seq_records:
                locus = seq_record.name if seq_record.name else seq_record.id
                feature_container = {}
                feature_counter = {}

                for feature in seq_record.features:
                    feature_type = feature.type
                    if feature_type not in feature_container:
                        feature_container[feature_type] = {}
                        feature_counter[feature_type] = 0

                    current_idx = feature_counter[feature_type]
                    feature_counter[feature_type] += 1

                    # Support simple and compound feature locations uniformly.
                    parts = feature.location.parts if hasattr(feature.location, "parts") else [feature.location]
                    locations = []
                    for part in parts:
                        try:
                            locations.append([int(part.start), int(part.end)])
                        except (TypeError, ValueError):
                            continue

                    strand = '-' if feature.location.strand == -1 else '+'
                    parsed_feature = {
                        'location': locations,
                        'strand': strand,
                    }

                    # Keep qualifier keys unchanged and flatten values to strings.
                    for qualifier_type, qualifier_values in feature.qualifiers.items():
                        if not qualifier_values:
                            continue
                        if isinstance(qualifier_values, list):
                            parsed_feature[qualifier_type] = qualifier_values[0] if len(qualifier_values) == 1 else ' '.join(str(x) for x in qualifier_values)
                        else:
                            parsed_feature[qualifier_type] = str(qualifier_values)

                    feature_container[feature_type][current_idx] = parsed_feature

                records[locus] = {
                    'locus': locus,
                    'features': feature_container,
                }

            return records

        def parse_gff(file_path) -> dict:
            """
            Parse a GFF3 (General Feature Format) file into a dictionary structure.

            :param file_path: path to genebank file
            :return: nested dictionary

            """
            records = {}
            with _get_line_iterator(file_path) as file:
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
            with _get_line_iterator(file_path) as file:
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
            annotation_type = detect_annotation_type(annotation)
        except ValueError as err:
            raise err

        # read in the annotation
        annotations = parse_functions[annotation_type](annotation)

        # sanity check whether one of the annotation ids and alignment ids match
        annotation_found = False
        for annotation in annotations.keys():
            for aln_id in aln:
                aln_id_sanitized = aln_id.split(' ')[0]
                # check in both directions
                if aln_id_sanitized in annotation:
                    annotation_found = True
                    break
                if annotation in aln_id_sanitized:
                    annotation_found = True
                    break

        if not annotation_found:
            raise ValueError(f'the annotations of {annotation} do not match any ids in the MSA')

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
