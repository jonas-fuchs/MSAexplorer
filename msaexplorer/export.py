"""
# Export module

This module lets you export data produced with MSA explorer.

## Functions:
"""

import os

from numpy import ndarray

from msaexplorer import config


def _check_and_create_path(path: str):
    """
    Check and create path if it doesn't exist.
    :param path: string to file
    """
    if path is not None:
        output_dir = os.path.dirname(path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)


def snps(snp_dict: dict, path: str | None = None, format_type: str = 'vcf') -> str | None | ValueError:
    """
    Export a SNP dictionary to a VCF or tabular format. Importantly, the input dictionary has to be in the standard
    format that MSAexplorer produces.

    :param snp_dict: Dictionary containing SNP positions and variant information.
    :param path: Path to output VCF or tabular format. (optional)
    :param format_type: Format type ('vcf' or 'tabular'). Default is 'vcf'.
    :return: A string containing the SNP data in the requested format.
    :raises ValueError: if the input dictionary is missing required keys or format_type is invalid.
    """

    def _validate():
        if not isinstance(snp_dict, dict):
            raise ValueError('Input SNP data must be a dictionary.')
        for key in ['#CHROM', 'POS']:
            if key not in snp_dict:
                raise ValueError(f"Missing required key '{key}' in SNP data.")
        if not isinstance(snp_dict['POS'], dict):
            raise ValueError('Expected the \'POS\' key to contain a dictionary of positions.')
        if format_type not in ['vcf', 'tabular']:
            raise ValueError('Invalid format_type.')
        _check_and_create_path(path)

    def _vcf_format(snp_dict: dict) -> list:
        """
        Produce  vcf formatted SNP data.
        :param snp_dict: dictionary containing SNP positions and variant information.
        :return: list of lines to write
        """
        output_lines = []
        # VCF header
        output_lines.append('##fileformat=VCFv4.2')
        output_lines.append('##source=MSAexplorer')
        output_lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
        # process each SNP position in sorted order
        for pos in sorted(snp_dict['POS'].keys()):
            pos_info = snp_dict['POS'][pos]
            ref = pos_info.get('ref', '.')
            alt_dict = pos_info.get('ALT', {})
            # Create comma-separated list of alternative alleles
            alt_alleles = ",".join(alt_dict.keys()) if alt_dict else "."
            # Prepare INFO field: include allele frequencies and sequence IDs
            afs = []
            seq_ids = []
            for alt, details in alt_dict.items():
                af = details.get('AF', 0)
                afs.append(str(af))
                seq_ids.append("|".join(details.get('SEQ_ID', [])))
            info_fields = []
            if afs:
                info_fields.append("AF=" + ",".join(afs))
            if seq_ids:
                info_fields.append("SEQ_ID=" + ",".join(seq_ids))
            info = ";".join(info_fields) if info_fields else "."

            # VCF is 1-indexed; we assume pos is 0-indexed and add 1
            line = f"{snp_dict['#CHROM']}\t{pos + 1}\t.\t{ref}\t{alt_alleles}\t.\t.\t{info}"
            output_lines.append(line)

        return output_lines

    def _tabular_format(snp_dict: dict) -> list:
        """
        Produce  tabular formatted SNP data.

        :param snp_dict: dictionary containing SNP positions and variant information.
        :return: list of lines to write
        """
        output_lines = []
        # Create a header for the tabular output
        output_lines.append('CHROM\tPOS\tREF\tALT\tAF\tSEQ_ID')

        # Process each SNP position and each alternative allele
        for pos in sorted(snp_dict['POS'].keys()):
            pos_info = snp_dict['POS'][pos]
            ref = pos_info.get('ref', '.')
            alt_dict = pos_info.get('ALT', {})
            for alt, details in alt_dict.items():
                af = details.get('AF', 0)
                seq_id = ",".join(details.get('SEQ_ID', []))
                output_lines.append(f"{snp_dict['#CHROM']}\t{pos + 1}\t{ref}\t{alt}\t{af}\t{seq_id}")

        return output_lines

    # validate correct input format
    _validate()

    # generate line data
    if format_type == 'vcf':
        lines = _vcf_format(snp_dict)
    else:
        lines = _tabular_format(snp_dict)

    # export to file or return plain text
    if path is not None:
        out_path = f"{path}.{format_type}"
        with open(out_path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)


def fasta(sequence: str, path: str | None = None,  header: str = 'consensus') -> str | None:
    """
    Export a fasta file to either a string or save directly to file.
    :param sequence: sequence to export
    :param path: path to save the file
    :param header: optional header file
    :return: fasta formatted string
    """
    def _validate_sequence(sequence: str):
        if not set(sequence).issubset(set(config.POSSIBLE_CHARS)):
            raise ValueError(f'Sequence contains invalid characters{set(sequence)}')

    _check_and_create_path(path)
    _validate_sequence(sequence)
    fasta_formated_sequence = f'>{header}\n{sequence}'
    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write(fasta_formated_sequence)
    else:
        return fasta_formated_sequence


def stats(stat_data: list | ndarray, seperator: str, path: str | None = None) -> str | None:
    """
    Export a list of stats per nucleotide to tabular or csv format.

    :param stat_data: list of stat values
    :param seperator: seperator for values and index
    :param path: path to save the file
    :return: tabular/csv formatted string
    """
    # ini
    _check_and_create_path(path)

    lines = [f'position{seperator}value']

    for idx, stat_val in enumerate(stat_data):
        lines.append(f'{idx}{seperator}{stat_val}')

    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)


def orf(orf_dict: dict, chrom: str, path: str | None = None) -> str | ValueError:
    """
    Exports the ORF dictionary to a .bed file.

    :param orf_dict: Dictionary containing ORF information.
    :param path: Path to the output .bed file.
    :param : Reference name
    """
    if not orf_dict:
        raise ValueError("The ORF dictionary is empty. Nothing to export.")
    else:
        if list(orf_dict[list(orf_dict.keys())[0]].keys()) != ['location', 'frame', 'strand', 'conservation', 'internal']:
            raise ValueError("The ORF dictionary has not the right format.")

    _check_and_create_path(path)

    lines = []

    for orf_id, orf_data in orf_dict.items():
        lines.append(
            f"{chrom}\t{orf_data['location'][0][0]}\t{orf_data['location'][0][1]}\t{orf_id}\t{orf_data['conservation']:.2f}\t{orf_data['strand']}"
        )

    if path is not None:
        with open(path, 'w') as out_file:
            out_file.write('\n'.join(lines))
    else:
        return '\n'.join(lines)
