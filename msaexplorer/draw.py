"""
contains the functions for drawing graphs
"""

# built-in
from typing import Callable, Dict

# MSAexplorer
from msaexplorer import explore
from msaexplorer import config

# libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like


def _validate_input_parameters(aln: explore.MSA, ax: plt.Axes):
    if not isinstance(aln, explore.MSA):
        raise ValueError('alignment has to be an MSA class. use explore.MSA() to read in alignment')
    if not isinstance(ax, plt.Axes):
        raise ValueError('ax has to be an matplotlib axis')


def identity_plot(aln: explore.MSA, ax: plt.Axes, show_seq_names: bool = False, custom_seq_names: tuple = (), reference_color = 'lightsteelblue', aln_colors:dict = config.ALN_COLORS, show_mask:bool = True, show_gaps:bool = True, fancy_gaps:bool = False, show_mismatches: bool = True, show_ambiguties: bool = False, show_x_label: bool = True, show_legend: bool = False):
    """
    generates an alignment overview plot
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param show_seq_names: if True, shows sequence names
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param aln_colors: dictionary containing colors: dict(0: color0, 1: color0, 2: color0, 3: color0) -> 0: match, 1: mismatch, 2: mask (N|X), 3: gap
    :param show_mask: whether to show N or X chars otherwise it will be shown as match or mismatch
    :param show_gaps: whether to show gaps otherwise it will be shown as match or mismatch
    :param fancy_gaps: show gaps with a small black bar
    :param show_mismatches: whether to show mismatches otherwise it will be shown as match
    :param show_ambiguties: whether to show non-N ambiguities -> only relevant for RNA/DNA sequences
    :param show_x_label: whether to show x label
    """

    def find_stretches(arr, target_value):
        """
        finds consecutive stretches of a number in an array
        :param arr: values in array
        :param target_value: target value to search for stretches
        :return: list of stretches (start value + length)
        """
        # Create a boolean array
        if np.isnan(target_value):
            is_target = np.isnan(arr)
        else:
            is_target = arr == target_value

        # Identify where changes occur in the boolean array
        change_points = np.where(np.diff(is_target) != 0)[0] + 1
        start_indices = np.insert(change_points, 0, 0)
        end_indices = np.append(change_points, len(arr))

        # Collect stretches where the boolean array is True
        stretches = [
            (int(start), int(end - start))
            for start, end in zip(start_indices, end_indices)
            if np.all(is_target[start:end])
        ]

        return stretches

    # input check
    _validate_input_parameters(aln, ax)

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError('reference color is not a color')
    if aln_colors.keys() != config.ALN_COLORS.keys():
        raise ValueError('configure your dictionary with 0 to 3 key and associated colors. See config.ALN_COLORS')
    for key in aln_colors.keys():
        for key2 in aln_colors[key]:
            if key2 not in ['type', 'color']:
                raise ValueError('configure your dictionary - see config.ALN_COLORS')
    if not all([is_color_like(aln_colors[x]['color']) for x in aln_colors.keys()]):
        raise ValueError('configure your dictionary with 0 to 3 key and associated colors. See config.ALN_COLORS')

    # validate custom names and set show names to True
    if custom_seq_names:
        show_seq_names = True
        if not isinstance(custom_seq_names, tuple):
            raise ValueError('configure your custom names list: custom_names=(name1, name2...)')
        if len(custom_seq_names) != len(aln.alignment.keys()):
            raise ValueError('length of sequences not equal to number of custom names')

    if fancy_gaps:
        show_gaps = True

    # get alignment and identity array
    identity_aln = aln.calc_identity_alignment(encode_mask=show_mask, encode_gaps=show_gaps, encode_mismatches=show_mismatches, encode_ambiguities=show_ambiguties)
    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom
    # define the y position of the first sequence
    y_position = len(aln.alignment) - 0.8
    detected_identity_values = {0}
    # ini collection for patches
    col = []
    for sequence, seq_name in zip(identity_aln, aln.alignment.keys()):
        # ini identity patches
        col.append(patches.Rectangle((zoom[0], y_position), zoom[1] - zoom[0],0.8,
                                               facecolor=reference_color if seq_name == aln.reference_id else aln_colors[0]['color']
                                               )
                             )
        for identity_value in [3, 2, 4, 1]: # first plot gaps, then mask, then ambiguities, then mismatches
            stretches = find_stretches(sequence, identity_value)
            if not stretches:
                continue
            # add values for legend
            detected_identity_values.add(identity_value)
            for stretch in stretches:
                col.append(
                    patches.Rectangle(
                        (stretch[0] + zoom[0], y_position),
                        stretch[1],
                        0.8,
                        color=aln_colors[identity_value]['color'],
                        linewidth=None
                    )
                )
                if not all([identity_value == 3, fancy_gaps]):
                    continue
                col.append(
                    patches.Rectangle(
                        (stretch[0] + zoom[0], y_position+0.375),
                        stretch[1],
                        0.05,
                        color='black',
                        linewidth=None
                    )
                )
        y_position -= 1

    # custom legend
    if show_legend:
        custom_legend = [ax.add_line(plt.Line2D([], [], color=aln_colors[val]['color'], marker='s' ,markeredgecolor='grey', linestyle='', markersize=10)) for val in
                         detected_identity_values]
        ax.legend(
            custom_legend,
            [aln_colors[x]['type'] for x in detected_identity_values],
            loc='upper right',
            bbox_to_anchor=(1, 1.1),
            ncols=len(detected_identity_values),
            frameon=False
        )

    # configure axis
    ax.add_collection(PatchCollection(col, match_original=True, linewidths='none'))
    ax.set_ylim(0, len(aln.alignment)+0.8/4)
    ax.set_xlim(zoom[0] - aln.length / 50, zoom[0] + aln.length + aln.length / 50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    if show_seq_names:
        ax.set_yticks(np.arange(len(aln.alignment)) + 0.6)
        if custom_seq_names:
            ax.set_yticklabels(custom_seq_names[::-1])
        else:
            ax.set_yticklabels(list(aln.alignment.keys())[::-1])
    else:
        ax.set_yticks([])
    if show_x_label:
        ax.set_xlabel('alignment position')


def stat_plot(aln: explore.MSA, ax: plt.Axes, stat_type: str, line_color: str = 'burlywood', rolling_average: int = 20, show_x_label: bool = False):
    """
    generate a plot for the various alignment stats
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param stat_type: 'entropy', 'gc' or 'coverage'
    :param line_color: color of the line
    :param rolling_average: average rolling window size in nucleotides or amino acids
    :param show_x_label: whether to show the x-axis label
    """

    def moving_average(arr, window_size):
        if window_size > 1:
            i = 0
            moving_averages, plotting_idx = [], []
            while i < len(arr) - window_size + 1:
                window = arr[i: i + window_size]
                moving_averages.append(sum(window) / window_size)
                plotting_idx.append(i + window_size / 2)
                i += 1

            return np.array(moving_averages), np.array(plotting_idx) if aln.zoom is None else np.array(plotting_idx) + aln.zoom[0]
        else:
            return arr, np.arange(aln.zoom[0], aln.length) if aln.zoom is not None else np.arange(aln.length)

    # define possible functions to calc here
    stat_functions: Dict[str, Callable[[], list]] = {
        'gc': aln.calc_gc,
        'entropy': aln.calc_entropy,
        'coverage': aln.calc_coverage
    }

    if stat_type not in stat_functions:
        raise ValueError('stat_type must be one of {}'.format(list(stat_functions.keys())))

    # input check
    _validate_input_parameters(aln, ax)
    if not is_color_like(line_color):
        raise ValueError('line color is not a color')

    # validate rolling average
    if rolling_average < 0 or rolling_average > aln.length:
        raise ValueError('rolling_average must be between 0 and length of sequence')

    # generate input data
    data, plot_idx = moving_average(stat_functions[stat_type](), rolling_average)
    # plot
    ax.plot(plot_idx, data, color=line_color)

    # specific visual cues for individual plots
    if stat_type == 'entropy' or stat_type == 'coverage':
        ax.fill_between(plot_idx, data, color=(line_color, 0.5))
    if stat_type == 'gc':
        ax.hlines(0.5, xmin=0, xmax=aln.zoom[0] + aln.length if aln.zoom is not None else aln.length, color='black', linestyles='-', linewidth=1)

    # format axis
    ax.set_ylim(0, 1.1)
    ax.set_xlim(
        (aln.zoom[0] - aln.length / 50, aln.zoom[0] + aln.length + aln.length / 50)
        if aln.zoom is not None else
        (0- aln.length / 50, aln.length + aln.length / 50)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if show_x_label:
        ax.set_xlabel('alignment position')
    if rolling_average > 1:
        ax.set_ylabel(f'{stat_type}\n rolling average')
    else:
        ax.set_ylabel(f'{stat_type}')


def variant_plot(aln: explore.MSA, ax: plt.Axes, show_x_label: bool = False, colors: dict|None = None, show_legend: bool = True):
    """
    Plots variants
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param show_x_label:  whether to show the x-axis label
    :param colors: colors for variants - standard are config.AA_colors or config.NT_colors
    :param legend: whether to show the legend or not
    """

    # validate input
    _validate_input_parameters(aln, ax)
    # define colors
    if colors is None:
        # define colors to use
        if aln.aln_type == 'AA':
            colors = config.AA_COLORS
        else:
            colors = config.NT_COLORS
    else:
        if colors is not dict:
            raise TypeError('Format colors like in config.AA_COLORS or config.NT_COLORS')
        if aln.aln_type == 'AA':
            if colors.keys() != config.AA_COLORS.keys():
                raise TypeError('Format colors like in config.AA_COLORS')
        if aln.aln_type == 'NT':
            if colors.keys() != config.NT_COLORS.keys():
                raise TypeError('Format colors like in config.NT_COLORS')
        for color in colors:
            if not is_color_like(color):
                raise TypeError(f'{color} is not a color')

    # get snps
    snps = aln.get_snps()
    # define where to plot (each ref type gets a separate line)
    ref_y_positions, y_pos, detected_var = {}, 0, set()

    # iterate over snp dict
    for pos in snps['POS']:
        for identifier in snps['POS'][pos]:
            # fill in y pos dict
            if identifier == 'ref':
                if snps['POS'][pos]['ref'] not in ref_y_positions:
                    ref_y_positions[snps['POS'][pos]['ref']] = y_pos
                    y_pos += 1
                    continue
            # plot
            if identifier == 'ALT':
                for alt in snps['POS'][pos]['ALT']:
                    ax.vlines(x=pos + aln.zoom[0] if aln.zoom is not None else pos,
                              ymin=ref_y_positions[snps['POS'][pos]['ref']],
                              ymax=ref_y_positions[snps['POS'][pos]['ref']] + snps['POS'][pos]['ALT'][alt]['AF'],
                              color=colors[alt],
                              zorder=100
                              )
                    ax.plot(pos + aln.zoom[0] if aln.zoom is not None else pos,
                            ref_y_positions[snps['POS'][pos]['ref']] + snps['POS'][pos]['ALT'][alt]['AF'],
                            color=colors[alt],
                            marker='o',
                            markersize=3)
                    detected_var.add(alt)
    # plot hlines
    for y_char in ref_y_positions:
        ax.hlines(
            ref_y_positions[y_char],
            xmin=aln.zoom[0] - aln.length / 50 if aln.zoom is not None else 0 - aln.length / 50,
            xmax=aln.zoom[0] + aln.length + aln.length / 50 if aln.zoom is not None else aln.length + aln.length / 50,
            color='black',
            linestyle='-',
            zorder=0,
            linewidth=0.75
        )
    # create a custom legend
    if show_legend:
        custom_legend = [ax.add_line(plt.Line2D([], [], color=colors[char], marker='o', linestyle='', markersize=5)) for char in
                         detected_var]
        ax.legend(
            custom_legend,
            detected_var,
            loc='upper right',
            title='variant',
            bbox_to_anchor=(1, 1.1),
            ncols=len(detected_var),
            frameon=False
        )

    # format axis
    ax.set_xlim(
        (aln.zoom[0] - aln.length / 50, aln.zoom[0] + aln.length + aln.length / 50)
        if aln.zoom is not None else
        (0 - aln.length / 50, aln.length + aln.length / 50)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([ref_y_positions[x] for x in ref_y_positions])
    ax.set_yticklabels(ref_y_positions.keys())
    ax.set_ylim(0, y_pos)
    if show_x_label:
        ax.set_xlabel('alignment position')
    ax.set_ylabel('reference')