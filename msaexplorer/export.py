"""
# Export module

This module lets you export data produced with MSA explorer.

## Functions:
"""

import numpy as np
from numpy import ndarray
from msaexplorer import config
from msaexplorer._data_classes import AlignmentStats, OrfCollection, VariantCollection
from msaexplorer._helpers import _check_and_create_path


def snps(snp_data: VariantCollection, format_type: str = 'vcf', path: str | None = None) -> str | None:
    """
    Export SNP data from a VariantCollection to VCF or tabular format.

    :param snp_data: VariantCollection containing SNP positions and variant information.
    :param format_type: Format type ('vcf' or 'tabular'). Default is 'vcf'.
    :param path: Path to output VCF or tabular format. (optional)
    :return: A string containing the SNP data in the requested format.
    :raises ValueError: If the input type is invalid or format_type is invalid.
    """

    def _validate():
        if not isinstance(snp_data, VariantCollection):
            raise ValueError('Input SNP data must be a VariantCollection dataclass.')
        if format_type not in ['vcf', 'tabular']:
            raise ValueError('Invalid format_type.')
        _check_and_create_path(path)

    def _vcf_format(data: VariantCollection) -> list[str]:
        """Produce VCF formatted SNP data."""
        output_lines = [
            '##fileformat=VCFv4.2',
            '##source=MSAexplorer',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
        ]

        for pos in sorted(data.positions.keys()):
            pos_info = data.positions[pos]
            alt_dict = pos_info.alt
            alt_alleles = ','.join(alt_dict.keys()) if alt_dict else '.'

            afs = [str(af) for af, _seq_ids in alt_dict.values()]
            seq_ids = ['|'.join(seq_ids) for _af, seq_ids in alt_dict.values()]
            info_fields = []
            if afs:
                info_fields.append('AF=' + ','.join(afs))
            if seq_ids:
                info_fields.append('SEQ_ID=' + ','.join(seq_ids))
            info = ';'.join(info_fields) if info_fields else '.'

            output_lines.append(
                f'{data.chrom}\t{pos + 1}\t.\t{pos_info.ref}\t{alt_alleles}\t.\t.\t{info}'
            )

        return output_lines

    def _tabular_format(data: VariantCollection) -> list[str]:
        """Produce tabular formatted SNP data."""
        output_lines = ['CHROM\tPOS\tREF\tALT\tAF\tSEQ_ID']

        for pos in sorted(data.positions.keys()):
            pos_info = data.positions[pos]
            for alt, (af, seq_ids) in pos_info.alt.items():
                output_lines.append(
                    f'{data.chrom}\t{pos + 1}\t{pos_info.ref}\t{alt}\t{af}\t{",".join(seq_ids)}'
                )

        return output_lines

    _validate()
    lines = _vcf_format(snp_data) if format_type == 'vcf' else _tabular_format(snp_data)

    if path is not None:
        out_path = f'{path}.{format_type}'
        with open(out_path, 'w') as out_file:
            out_file.write('\n'.join(lines))
        return None

    return '\n'.join(lines)


def fasta(sequence: str | dict, header: str | None = None, path: str | None = None) -> str | None:
    """
    Export a fasta sequence from str or alignment in dictionary format to either a string or save directly to file.
    The alignment format must have headers as keys and the corresponding sequence as values.
    :param sequence: sequence to export
    :param header: optional header file
    :param path: path to save the file
    :return: fasta formatted string
    """
    def _validate_sequence(seq: str):
        if not set(seq).issubset(set(config.POSSIBLE_CHARS)):
            raise ValueError(f'Sequence contains invalid characters. Detected chars: {set(seq)}')

    _check_and_create_path(path)
    fasta_formated_sequence = ''

    if type(sequence) is str:
        _validate_sequence(sequence)
        fasta_formated_sequence = f'>{header}\n{sequence}'
    elif type(sequence) is dict:
        for header, sequence in sequence.items():
            if type(sequence) is not str:
                raise ValueError(f'Sequences in the dictionary must be strings.')
            _validate_sequence(sequence)
            fasta_formated_sequence = f'{fasta_formated_sequence}\n>{header}\n{sequence}' if fasta_formated_sequence != '' else f'>{header}\n{sequence}'

    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write(fasta_formated_sequence)
    else:
        return fasta_formated_sequence


def stats(stat_data: AlignmentStats | list | ndarray, seperator: str = '\t', path: str | None = None) -> str | None:
    """
    Export a list of stats per nucleotide to tabular or csv format.

    :param stat_data: position statistic dataclass or list/array of values
    :param seperator: seperator for values and index
    :param path: path to save the file
    :return: tabular/csv formatted string
    """
    # ini
    _check_and_create_path(path)

    lines = [f'position{seperator}value']

    if isinstance(stat_data, AlignmentStats):
        positions = stat_data.positions
        values = stat_data.values
    else:
        values = stat_data
        positions = np.arange(len(values), dtype=int)

    for position, stat_val in zip(positions, values):
        lines.append(f'{position}{seperator}{stat_val}')

    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)


def orf(orfs: OrfCollection, chrom: str, path: str | None = None) -> str | ValueError | None:
    """
    Exports the ORF collection to a .bed file.

    :param orf_dict: OrfContainer instance
    :param chrom: CHROM identifier for bed format.
    :param path: Path to the output .bed file.
    """
    if not isinstance(orfs, OrfCollection):
        raise ValueError('The ORF collection must be an instance of msaexplorer._data_classes.OrfCollection.')

    if not orfs:
        raise ValueError('The ORF collection is empty.')

    _check_and_create_path(path)

    lines = []

    for orf_id, orf_data in orfs.items():
        loc = orf_data.location[0]
        conservation = orf_data.conservation
        strand = orf_data.strand
        lines.append(f"{chrom}\t{loc[0]}\t{loc[1]}\t{orf_id}\t{conservation:.2f}\t{strand}")

    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
        return None
    else:
        return '\n'.join(lines)


def character_freq(char_dict: dict, seperator: str = '\t', path: str | None = None) -> str | None | ValueError:
    """
    Export a character frequency dictionary to tabular or csv format.

    :param char_dict: Dictionary containing the character frequencies.
    :param seperator: seperator for the table e.g. tabular or comma
    :param path: Path to output table.

    :return: A string containing the character frequency table.
    :raises ValueError: if the input dictionary is missing required keys or format_type is invalid.
    """

    def _validate():
        if not isinstance(char_dict, dict):
            raise ValueError('Data must be a dictionary.')
        for key, value in char_dict.items():
            for key_2, value_2 in value.items():
                if key_2 not in config.POSSIBLE_CHARS:
                    raise ValueError(f'The key {key_2} is invalid.')
                for key_3, value_3 in value_2.items():
                    if key_3 not in ['counts', '% of alignment', '% of non-gapped']:
                        raise ValueError(f'The key "{key_3}" is invalid.')

    # validate input
    _validate()

    lines = [F'sequence{seperator}char{seperator}counts{seperator}% of non-gapped']
    for key, value in char_dict.items():
        if key == 'total':
            continue
        for key_2, value_2 in value.items():
            if key_2 == '-':
                continue
            lines.append(f'{key}{seperator}{key_2}{seperator}{value_2["counts"]}{seperator}{value_2["% of non-gapped"]}')

    # export data
    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)


def percent_recovery(rec_dict: dict, seperator: str = '\t', path: str | None = None) -> str | None | ValueError:
    """
    Export percent_recovery dictionary to tabular or csv format.

    :param rec_dict: Dictionary containing the character frequencies.
    :param seperator: seperator for the table e.g. tabular or comma
    :param path: Path to output table.

    :return: A string containing the character frequency table.
    :raises ValueError: if the input dictionary is missing required keys or format_type is invalid.
    """
    def _validate():
        if not isinstance(rec_dict, dict):
            raise ValueError('Data must be a dictionary.')
        for key, value in rec_dict.items():
            if type(key) != str:
                raise ValueError(f'The key {key} is invalid.')
            elif type(value) != float:
                raise ValueError(f'The value {value} is invalid.')

    # validate input
    _validate()

    lines = [F'sequence{seperator}% recovery']
    for key, value in rec_dict.items():
        lines.append(
            f'{key}{seperator}{value}')

    # export data
    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)

