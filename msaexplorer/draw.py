"""
contains the functions for drawing graphs
"""

# built-in
from typing import Callable, Dict
from copy import deepcopy

# MSAexplorer
from msaexplorer import explore
from msaexplorer import config

# libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.cm import ScalarMappable
from matplotlib.colors import is_color_like, Normalize
from matplotlib.collections import PatchCollection


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
    ax.set_xlim(zoom[0], zoom[0] + aln.length)
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
    :param rolling_average: average rolling window size left and right of a position in nucleotides or amino acids
    :param show_x_label: whether to show the x-axis label
    """

    def moving_average(arr, window_size):
        if window_size > 1:
            i = 0
            moving_averages, plotting_idx = [], []
            while i < len(arr) + 1:
                window_left = arr[i- window_size: i] if i > window_size else arr[0:i]
                window_right = arr[i: i + window_size] if i < len(arr) - window_size else arr[i:len(arr)]
                moving_averages.append((sum(window_left)+ sum(window_right)) / (len(window_left) + len(window_right)))
                plotting_idx.append(i)
                i += 1

            return np.array(moving_averages), np.array(plotting_idx) if aln.zoom is None else np.array(plotting_idx) + aln.zoom[0]
        else:
            return arr, np.arange(aln.zoom[0], aln.zoom[1]) if aln.zoom is not None else np.arange(aln.length)

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
        ax.hlines(0.5, xmin=0, xmax=aln.zoom[0] + aln.length if aln.zoom is not None else aln.length, color='black', linestyles='--', linewidth=1)

    # format axis
    ax.set_ylim(0, 1.1)
    ax.set_xlim(
        (aln.zoom[0], aln.zoom[0] + aln.length)
        if aln.zoom is not None else
        (0, aln.length)
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
    :param show_legend: whether to show the legend
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
            xmin=aln.zoom[0] if aln.zoom is not None else 0,
            xmax=aln.zoom[0] + aln.length if aln.zoom is not None else aln.length,
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
            bbox_to_anchor=(1, 1.15),
            ncols=len(detected_var),
            frameon=False
        )

    # format axis
    ax.set_xlim(
        (aln.zoom[0], aln.zoom[0] + aln.length)
        if aln.zoom is not None else
        (0, aln.length)
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


def orf_plot(aln: explore.MSA, ax: plt.Axes, min_length: int = 500, non_overlapping_orfs: bool = True, cmap: str = 'Blues', show_direction:bool = True, show_x_label: bool = False, show_legend: bool = True):
    """
    Plot conserved ORFs
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param min_length: minimum length of orf
    :param non_overlapping_orfs: whether to consider overlapping orfs
    :param cmap: color mapping for % identity - see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    :param show_direction: show strand information
    :param show_x_label: whether to show the x-axis label
    :param show_legend: whether to show the legend
    """
    # helper function
    def add_track_positions(annotation_dic):
        # create a dict and sort
        annotation_dic = dict(sorted(annotation_dic.items(), key=lambda x: x[1]['positions'][0]))

        # remember for each track the largest stop
        track_stops = [0]

        for ann in annotation_dic:
            track = 0
            # check if a start of a gene is smaller than the stop of the current track
            # -> move to new track
            while annotation_dic[ann]['positions'][0] < track_stops[track]:
                track += 1
                # if all prior tracks are potentially causing an overlap
                # create a new track and break
                if len(track_stops) <= track:
                    track_stops.append(0)
                    break
            # in the current track remember the stop of the current gene
            track_stops[track] = annotation_dic[ann]['positions'][1]
            # and indicate the track in the dict
            annotation_dic[ann]['track'] = track

        return annotation_dic

    cmap = ScalarMappable(norm=Normalize(0, 100), cmap=plt.get_cmap('Blues'))

    # validate input
    _validate_input_parameters(aln, ax)

    # get orfs --> first deepcopy and reset zoom that the orfs are also zoomed in (by default, the orfs are only
    # searched within the zoomed region)
    aln_temp = deepcopy(aln)
    aln_temp.zoom = None
    if non_overlapping_orfs:
        annotation_dict = add_track_positions(aln_temp.get_non_overlapping_conserved_orfs(min_length=min_length))
    else:
        annotation_dict = add_track_positions(aln_temp.get_conserved_orfs(min_length=min_length))
    # plot
    for annotation in annotation_dict:
        x_value = annotation_dict[annotation]['positions'][0] + aln.zoom[0] if aln.zoom is not None else annotation_dict[annotation]['positions'][0]
        length = annotation_dict[annotation]['positions'][1] - annotation_dict[annotation]['positions'][0]
        ax.add_patch(
            patches.FancyBboxPatch(
                (x_value, annotation_dict[annotation]['track'] + 0.2),
                length,
                0.8,
                boxstyle="round,pad=-0.0040,rounding_size=0.1",
                ec="black",
                fc=cmap.to_rgba(annotation_dict[annotation]['conservation'])
            )
        )
        if show_direction:
            if annotation_dict[annotation]['strand'] == '-':
                marker = '<'
            else:
                marker = '>'
            ax.plot(x_value + length/2, annotation_dict[annotation]['track'] + 0.6, marker=marker, markersize=5, color='white', markeredgecolor='black')

    # legend
    if show_legend:
        cbar = plt.colorbar(cmap, ax=ax, location= 'top', orientation='horizontal', fraction=0.2, pad=0.1, anchor=(1,0))
        cbar.set_label('% identity')
        cbar.set_ticks([0, 100])

    # format axis
    ax.set_xlim(
        (aln.zoom[0], aln.zoom[0] + aln.length)
        if aln.zoom is not None else
        (0, aln.length)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    if show_x_label:
        ax.set_xlabel('alignment position')
    ax.set_ylim(bottom=0)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_title('conserved orfs', loc='left')

