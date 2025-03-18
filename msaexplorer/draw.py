r"""
# The draw module

The draw module lets you draw alignments and statistic plots such as SNPs, ORFs, entropy and much more. For each plot a
`matplotlib axes` has to be passed to the plotting function.

Importantly some of the plotting features can only be accessed for nucleotide alignments but not for amino acid alignments.
The functions will raise the appropriate exception in such a case.

## Functions

"""

# built-in
from itertools import chain
from typing import Callable, Dict
from copy import deepcopy
import os

import matplotlib
from numpy import ndarray

# MSAexplorer
from msaexplorer import explore, config

# libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.cm import ScalarMappable
from matplotlib.colors import is_color_like, Normalize, to_rgba
from matplotlib.collections import PatchCollection


# general helper functions
def _validate_input_parameters(aln: explore.MSA, ax: plt.Axes, annotation: explore.Annotation | None = None):
    """
    Validate MSA class and axis.
    """
    if not isinstance(aln, explore.MSA):
        raise ValueError('alignment has to be an MSA class. use explore.MSA() to read in alignment')
    if not isinstance(ax, plt.Axes):
        raise ValueError('ax has to be an matplotlib axis')
    if annotation is not None:
        if not isinstance(annotation, explore.Annotation):
            raise ValueError('annotation has to be an annotation class. use explore.Annotation() to read in annotation')


def _format_x_axis(aln: explore.MSA, ax: plt.Axes, show_x_label: bool, show_left: bool):
    """
    General axis formatting.
    """
    ax.set_xlim(
        (aln.zoom[0] - 0.5, aln.zoom[0] + aln.length - 0.5) if aln.zoom is not None else (-0.5, aln.length - 0.5)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if show_x_label:
        ax.set_xlabel('alignment position')
    if not show_left:
        ax.spines['left'].set_visible(False)


# helper functions for alignment plots
def _find_stretches(arr, target):
    """
    Finds consecutive stretches of a number or np.nan in an array.
    :param arr: values in array
    :param target: target value to search for stretches
    :return: list of stretches (start value + length)
    """

    stretches = []

    # Create a boolean array
    if np.isnan(target):
        matches = np.isnan(arr)
    else:
        matches = arr == target

    start_index = None  # To track the start of a stretch

    for i, is_match in enumerate(matches):
        if is_match:
            if start_index is None:  # Start of a new stretch
                start_index = i
        else:
            if start_index is not None:  # End of a stretch
                stretch_length = i - start_index
                stretches.append((start_index, stretch_length))
                start_index = None

    # If the last stretch reaches the end of the array
    if start_index is not None:
        stretch_length = len(matches) - start_index
        stretches.append((start_index, stretch_length))

    return stretches


def _seq_names(aln: explore.MSA, ax: plt.Axes, custom_seq_names: tuple, show_seq_names: bool):
    """
    Validate custom names and set show names to True. Format axis accordingly.
    """
    if custom_seq_names:
        show_seq_names = True
        if not isinstance(custom_seq_names, tuple):
            raise ValueError('configure your custom names list: custom_names=(name1, name2...)')
        if len(custom_seq_names) != len(aln.alignment.keys()):
            raise ValueError('length of sequences not equal to number of custom names')
    if show_seq_names:
        ax.yaxis.set_ticks_position('none')
        ax.set_yticks(np.arange(len(aln.alignment)) + 0.6)
        if custom_seq_names:
            ax.set_yticklabels(custom_seq_names[::-1])
        else:
            names = [x.split(' ')[0] for x in list(aln.alignment.keys())[::-1]]
            ax.set_yticklabels(names)
    else:
        ax.set_yticks([])


def _create_identity_patch(aln: explore.MSA, col: list, zoom: tuple[int, int], y_position: float | int, reference_color: str, seq_name: str, identity_color: str | ndarray):
    """
    Creates the initial patch.
    """

    col.append(patches.Rectangle((zoom[0] - 0.5, y_position), zoom[1] - zoom[0], 0.8,
                                                   facecolor=reference_color if seq_name == aln.reference_id else identity_color
                                                   )
                                 )


def _plot_annotation(annotation_dict: dict, ax: plt.Axes, direction_marker_size: int | None, color: str | ScalarMappable):
    """
    Plot annotations
    :param annotation_dict: dict of annotations
    :param ax: matplotlib Axes
    :param direction_marker_size: size of marker
    :param color: color of annotation (color or scalar)
    """
    for annotation in annotation_dict:
        for locations in annotation_dict[annotation]['location']:
            x_value = locations[0]
            length = locations[1] - locations[0]
            ax.add_patch(
                patches.FancyBboxPatch(
                    (x_value, annotation_dict[annotation]['track'] + 1),
                    length,
                    0.8,
                    boxstyle="Round, pad=0",
                    ec="black",
                    fc=color.to_rgba(annotation_dict[annotation]['conservation']) if isinstance(color, ScalarMappable) else color,
                )
            )
            if direction_marker_size is not None:
                if annotation_dict[annotation]['strand'] == '-':
                    marker = '<'
                else:
                    marker = '>'
                ax.plot(x_value + length/2, annotation_dict[annotation]['track'] + 1.4, marker=marker, markersize=direction_marker_size, color='white', markeredgecolor='black')

        # plot linked annotations (such as splicing)
        if len(annotation_dict[annotation]['location']) > 1:
            y_value = annotation_dict[annotation]['track'] + 1.4
            start = None
            for locations in annotation_dict[annotation]['location']:
                if start is None:
                    start = locations[1]
                    continue
                ax.plot([start, locations[0]], [y_value, y_value], '--', linewidth=2, color='black')
                start = locations[1]


def _create_stretch_patch(col: list, stretches: list, zoom: tuple[int, int], y_position: float | int, colors: dict | ndarray, fancy_gaps: bool, matrix_value: int | float | ndarray):
    """
    Create a patch and add to list.
    """
    for stretch in stretches:
        col.append(
            patches.Rectangle(
                (stretch[0] - 0.5 + zoom[0], y_position),
                stretch[1],
                0.8,
                color=colors,
                linewidth=None
            )
        )
        if not fancy_gaps:
            continue
        if not np.isnan(matrix_value):
            continue
        col.append(
            patches.Rectangle(
                (stretch[0] - 0.5 + zoom[0], y_position + 0.375),
                stretch[1],
                0.05,
                color='black',
                linewidth=None
            )
        )


def _add_track_positions(annotation_dic):
    # create a dict and sort
    annotation_dic = dict(sorted(annotation_dic.items(), key=lambda x: x[1]['location'][0][0]))

    # remember for each track the largest stop
    track_stops = [0]

    for ann in annotation_dic:
        flattened_locations = list(chain.from_iterable(annotation_dic[ann]['location']))  # flatten list
        track = 0
        # check if a start of a gene is smaller than the stop of the current track
        # -> move to new track
        while flattened_locations[0] < track_stops[track]:
            track += 1
            # if all prior tracks are potentially causing an overlap
            # create a new track and break
            if len(track_stops) <= track:
                track_stops.append(0)
                break
        # in the current track remember the stop of the current gene
        track_stops[track] = flattened_locations[-1]
        # and indicate the track in the dict
        annotation_dic[ann]['track'] = track

    return annotation_dic


def _get_contrast_text_color(rgba_color):
    """
    compute the brightness of a color
    """
    r, g, b, a = rgba_color
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return 'white' if brightness < 0.4 else 'black'


def _plot_sequence_text(aln, seq_name, ref_name, values, matrix, ax, zoom, y_position, value_to_skip, ref_color, cmap: None | ScalarMappable = None):

    x_text = 0
    if seq_name == ref_name:
        different_cols = np.any((matrix != value_to_skip) & ~np.isnan(matrix), axis=0)
    else:
        different_cols = [False]*aln.length

    for idx, (character, value) in enumerate(zip(aln.alignment[seq_name], values)):
        if value != value_to_skip and character != '-' or seq_name == ref_name and character != '-':
            # text color
            if seq_name == ref_name:
                text_color = _get_contrast_text_color(to_rgba(ref_color))
            elif cmap is not None:
                text_color = _get_contrast_text_color(cmap.to_rgba(value))
            else:
                text_color = 'black'

            ax.text(
                x=x_text + zoom[0] if zoom is not None else x_text,
                y=y_position + 0.4,
                s=character,
                fontweight='bold' if different_cols[idx] else 'normal',
                ha='center',
                va='center',
                c='grey' if not different_cols[idx] and seq_name == ref_name else text_color,
            )
        x_text += 1

def identity_alignment(aln: explore.MSA, ax: plt.Axes, show_title: bool = True, show_sequence: bool = False, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', show_mask:bool = True, show_gaps:bool = True, fancy_gaps:bool = False, show_mismatches: bool = True, show_ambiguities: bool = False, color_mismatching_chars: bool = False, show_x_label: bool = True, show_legend: bool = False, bbox_to_anchor: tuple[float|int, float|int] | list[float|int, float|int]= (1, 1)):
    """
    Generates an identity alignment overview plot.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param show_title: whether to show title
    :param show_seq_names: whether to show seq names
    :param show_sequence: whether to show sequence for differences and reference - zoom in to avoid plotting issues
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param show_mask: whether to show N or X chars otherwise it will be shown as match or mismatch
    :param show_gaps: whether to show gaps otherwise it will be shown as match or mismatch
    :param fancy_gaps: show gaps with a small black bar
    :param show_mismatches: whether to show mismatches otherwise it will be shown as match
    :param show_ambiguities: whether to show non-N ambiguities -> only relevant for RNA/DNA sequences
    :param color_mismatching_chars: color mismatching chars with their unique color
    :param show_x_label: whether to show x label
    :param show_legend: whether to show the legend
    :param bbox_to_anchor: bounding box coordinates for the legend - see: https://matplotlib.org/stable/api/legend_api.html
    """

    # input check
    _validate_input_parameters(aln, ax)

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError(f'{reference_color} for reference is not a color')

    # Both options for gaps work hand in hand
    if fancy_gaps:
        show_gaps = True

    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom

    # define colors and identity values (correspond to values in the alignment matrix)
    aln_colors = config.IDENTITY_COLORS
    if color_mismatching_chars:
        identity_values = [np.nan, -2, -3] if show_gaps else [-2, -3]
        colors_to_extend = config.CHAR_COLORS[aln.aln_type]
        identity_values = identity_values + [x+1 for x in list(range(len(colors_to_extend)))]  # x+1 is needed to allow correct mapping
        for idx, char in enumerate(colors_to_extend):
            aln_colors[idx + 1] = {'type': char, 'color': colors_to_extend[char]}
    else:
        identity_values = [np.nan, -1, -2, -3]

    # get alignment and identity array
    identity_aln = aln.calc_identity_alignment(
        encode_mask=show_mask,
        encode_gaps=show_gaps,
        encode_mismatches=show_mismatches,
        encode_ambiguities=show_ambiguities,
        encode_each_mismatch_char=color_mismatching_chars
    )

    # define the y position of the first sequence
    y_position = len(aln.alignment) - 0.8
    detected_identity_values = {0}

    # ini collection for patches
    col = []
    for values, seq_name in zip(identity_aln, aln.alignment.keys()):
        # ini identity patches
        _create_identity_patch(aln, col, zoom, y_position, reference_color, seq_name, aln_colors[0]['color'])
        # find and plot stretches that are different
        for identity_value in identity_values:
            stretches = _find_stretches(values, identity_value)
            if not stretches:
                continue
            # add values for legend
            detected_identity_values.add(identity_value)
            _create_stretch_patch(col, stretches, zoom, y_position, aln_colors[identity_value]['color'], fancy_gaps, identity_value)

        if show_sequence:
            _plot_sequence_text(aln, seq_name, aln.reference_id, values, identity_aln, ax, zoom, y_position, 0, reference_color)

        y_position -= 1

    # custom legend
    if show_legend:
        # create it
        custom_legend = [
            ax.add_line(
                plt.Line2D(
                    [],
                    [],
                    color=aln_colors[x]['color'], marker='s' ,markeredgecolor='grey', linestyle='', markersize=10)) for x in aln_colors if x in detected_identity_values
        ]
        # plot it
        ax.legend(
            custom_legend,
            [aln_colors[x]['type'] for x in aln_colors if x in detected_identity_values],
            loc='lower right',
            bbox_to_anchor=bbox_to_anchor,
            ncols=len(detected_identity_values) / 2 if aln.aln_type == 'AA' and color_mismatching_chars else len(detected_identity_values),
            frameon=False
        )

    # format seq names
    _seq_names(aln, ax, custom_seq_names, show_seq_names)

    # configure axis
    ax.add_collection(PatchCollection(col, match_original=True, linewidths='none', joinstyle='miter', capstyle='butt'))
    ax.set_ylim(0, len(aln.alignment))
    if show_title:
        ax.set_title('identity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def similarity_alignment(aln: explore.MSA, ax: plt.Axes, matrix_type: str | None = None, show_title: bool = True, show_sequence: bool = False, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', cmap: str = 'PuBu_r', gap_color: str = 'white', show_gaps:bool = True, fancy_gaps:bool = False, show_x_label: bool = True, show_cbar: bool = False, cbar_fraction: float = 0.1):
    """
    Generates a similarity alignment overview plot. Importantly the similarity values are normalized!
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param matrix_type: substitution matrix - see config.SUBS_MATRICES, standard: NT - TRANS, AA - BLOSUM65
    :param show_title: whether to show title
    :param show_sequence: whether to show sequence for differences and reference - zoom in to avoid plotting issues
    :param show_seq_names: whether to show seq names
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param cmap: color mapping for % identity - see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    :param gap_color: color for gaps
    :param show_gaps: whether to show gaps otherwise it will be ignored
    :param fancy_gaps: show gaps with a small black bar
    :param show_x_label: whether to show x label
    :param show_cbar: whether to show the legend - see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
    :param cbar_fraction: fraction of the original ax reserved for the legend
    """
    # input check
    _validate_input_parameters(aln, ax)

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError(f'{reference_color} for reference is not a color')

    # Both options for gaps work hand in hand
    if fancy_gaps:
        show_gaps = True

    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom

    # get data
    similarity_aln = aln.calc_similarity_alignment(matrix_type=matrix_type)  # use normalized values here
    similarity_aln = similarity_aln.round(2)  # round data for good color mapping

    # determine min max values of the underlying matrix and create cmap
    min_value, max_value = 0, 1
    cmap = ScalarMappable(
        norm=Normalize(
            vmin=min_value,
            vmax=max_value
        ),
        cmap=plt.get_cmap(cmap)
    )

    # create similarity values
    similarity_values = np.append(
        [np.nan], np.arange(start=min_value, stop=max_value, step=0.01), 0
    ) if show_gaps else np.arange(
        start=min_value, stop=max_value, step=0.01
    )
    # round it to be absolutely sure that values match with rounded sim alignment
    similarity_values = similarity_values.round(2)


    # create plot
    col = []
    y_position = len(aln.alignment) - 0.8

    for values, seq_name in zip(similarity_aln, aln.alignment.keys()):
        # create initial patch for each sequence
        _create_identity_patch(
            aln, col, zoom, y_position, reference_color, seq_name,
            cmap.to_rgba(max_value) if seq_name != aln.reference_id else reference_color
        )
        # find and plot stretches
        for similarity_value in similarity_values:
            # skip patch plotting for reference except if it is a gap
            if seq_name == aln.reference_id and not np.isnan(similarity_value):
                continue
            stretches = _find_stretches(values, similarity_value)
            # create stretches
            if stretches:
                _create_stretch_patch(
                    col,
                    stretches,
                    zoom,
                    y_position,
                    cmap.to_rgba(similarity_value) if not np.isnan(similarity_value) else gap_color,
                    fancy_gaps,
                    similarity_value
                )
        if show_sequence:
            _plot_sequence_text(aln, seq_name, aln.reference_id, values, similarity_aln, ax, zoom, y_position, 1, reference_color, cmap=cmap)

        # new y position for the next sequence
        y_position -= 1

    # legend
    if show_cbar:
        cbar = plt.colorbar(cmap, ax=ax, location= 'top', anchor=(1,0), shrink=0.2, pad=2/ax.bbox.height, fraction=cbar_fraction)
        cbar.set_ticks([min_value, max_value])
        cbar.set_ticklabels(['low', 'high'])

    # format seq names
    _seq_names(aln, ax, custom_seq_names, show_seq_names)

    # configure axis
    ax.add_collection(PatchCollection(col, match_original=True, linewidths='none'))
    ax.set_ylim(0, len(aln.alignment))
    if show_title:
        ax.set_title('similarity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def stat_plot(aln: explore.MSA, ax: plt.Axes, stat_type: str, line_color: str = 'burlywood', line_width: int | float = 2, rolling_average: int = 20, show_x_label: bool = False, show_title: bool = True):
    """
    Generate a plot for the various alignment stats.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param stat_type: 'entropy', 'gc', 'coverage', 'ts/tv', 'identity' or 'similarity' -> (here default matrices are used NT - TRANS, AA - BLOSUM65)
    :param line_color: color of the line
    :param line_width: width of the line
    :param rolling_average: average rolling window size left and right of a position in nucleotides or amino acids
    :param show_x_label: whether to show the x-axis label
    :param show_title: whether to show the title
    """

    def moving_average(arr, window_size):
        if window_size > 1:
            i = 0
            moving_averages, plotting_idx = [], []
            while i < len(arr) + 1:
                half_window_size = window_size // 2
                window_left = arr[i - half_window_size : i] if i > half_window_size else arr[0:i]
                window_right = arr[i: i + half_window_size] if i < len(arr) - half_window_size else arr[i: len(arr)]
                moving_averages.append((sum(window_left) + sum(window_right)) / (len(window_left) + len(window_right)))
                plotting_idx.append(i)
                i += 1

            return np.array(moving_averages), np.array(plotting_idx) if aln.zoom is None else np.array(plotting_idx) + aln.zoom[0]
        else:
            return arr, np.arange(aln.zoom[0], aln.zoom[1]) if aln.zoom is not None else np.arange(aln.length)

    # define possible functions to calc here
    stat_functions: Dict[str, Callable[[], list | ndarray]] = {
        'gc': aln.calc_gc,
        'entropy': aln.calc_entropy,
        'coverage': aln.calc_coverage,
        'identity': aln.calc_identity_alignment,
        'similarity': aln.calc_similarity_alignment,
        'ts tv score': aln.calc_transition_transversion_score
    }

    if stat_type not in stat_functions:
        raise ValueError('stat_type must be one of {}'.format(list(stat_functions.keys())))

    # input check
    _validate_input_parameters(aln, ax)
    if not is_color_like(line_color):
        raise ValueError('line color is not a color')

    # validate rolling average
    if rolling_average < 1 or rolling_average > aln.length:
        raise ValueError('rolling_average must be between 1 and length of sequence')

    # generate input data
    array = stat_functions[stat_type]()

    if stat_type == 'identity':
        min_value, max_value = -1, 0
    elif stat_type == 'ts tv score':
        min_value, max_value = -1, 1
    else:
        min_value, max_value = 0, 1
    if stat_type in ['identity', 'similarity']:
        # for the mean nan values get handled as the lowest possible number in the matrix
        array = np.nan_to_num(array, True, min_value)
        array = np.mean(array, axis=0)
    data, plot_idx = moving_average(array, rolling_average)

    # plot the data
    ax.fill_between(
        # this add dummy data left and right for better plotting
        # otherwise only half of the step is shown
        np.concatenate(([plot_idx[0] - 0.5], plot_idx, [plot_idx[-1] + 0.5])) if rolling_average == 1 else plot_idx,
        np.concatenate(([data[0]], data, [data[-1]])) if rolling_average == 1 else data,
        min_value,
        linewidth = line_width,
        edgecolor=line_color,
        step='mid' if rolling_average == 1 else None,
        facecolor=(line_color, 0.6) if stat_type not in ['ts tv score', 'gc'] else 'none'
    )
    if stat_type == 'gc':
        ax.hlines(0.5, xmin=0, xmax=aln.zoom[0] + aln.length if aln.zoom is not None else aln.length, color='black', linestyles='--', linewidth=1)

    # format axis
    ax.set_ylim(min_value, max_value*0.1+max_value)
    ax.set_yticks([min_value, max_value])
    if stat_type == 'gc':
        ax.set_yticklabels(['0', '100'])
    elif stat_type == 'ts tv score':
        ax.set_yticklabels(['tv', 'ts'])
    else:
        ax.set_yticklabels(['low', 'high'])

    # show title
    if show_title:
        ax.set_title(
            f'{stat_type} (average over {rolling_average} positions)' if rolling_average > 1 else f'{stat_type} for each position',
            loc='left'
        )

    _format_x_axis(aln, ax, show_x_label, show_left=True)


def variant_plot(aln: explore.MSA, ax: plt.Axes, lollisize: tuple[int, int] | list[int, int] = (1, 3), show_x_label: bool = False, show_legend: bool = True, bbox_to_anchor: tuple[float|int, float|int] | list[float|int, float|int] = (1, 1)):
    """
    Plots variants.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param lollisize: (stem_size, head_size)
    :param show_x_label:  whether to show the x-axis label
    :param show_legend: whether to show the legend
    :param bbox_to_anchor: bounding box coordinates for the legend - see: https://matplotlib.org/stable/api/legend_api.html
    """

    # validate input
    _validate_input_parameters(aln, ax)
    if not isinstance(lollisize, tuple) or len(lollisize) != 2:
        raise ValueError('lollisize must be tuple of length 2 (stem, head)')
    for _size in lollisize:
        if not isinstance(_size, float | int) or _size <= 0:
            raise ValueError('lollisize must be floats greater than zero')

    # define colors
    colors = config.CHAR_COLORS[aln.aln_type]
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
                    y_pos += 1.1
                    continue
            # plot
            if identifier == 'ALT':
                for alt in snps['POS'][pos]['ALT']:
                    ax.vlines(x=pos + aln.zoom[0] if aln.zoom is not None else pos,
                              ymin=ref_y_positions[snps['POS'][pos]['ref']],
                              ymax=ref_y_positions[snps['POS'][pos]['ref']] + snps['POS'][pos]['ALT'][alt]['AF'],
                              color=colors[alt],
                              zorder=100,
                              linewidth=lollisize[0]
                              )
                    ax.plot(pos + aln.zoom[0] if aln.zoom is not None else pos,
                            ref_y_positions[snps['POS'][pos]['ref']] + snps['POS'][pos]['ALT'][alt]['AF'],
                            color=colors[alt],
                            marker='o',
                            markersize=lollisize[1]
                            )
                    detected_var.add(alt)

    # plot hlines
    for y_char in ref_y_positions:
        ax.hlines(
            ref_y_positions[y_char],
            xmin=aln.zoom[0] - 0.5 if aln.zoom is not None else -0.5,
            xmax=aln.zoom[0] + aln.length + 0.5 if aln.zoom is not None else aln.length + 0.5,
            color='black',
            linestyle='-',
            zorder=0,
            linewidth=0.75
        )
    # create a custom legend
    if show_legend:
        custom_legend = [
            ax.add_line(
                plt.Line2D(
                    [],
                    [],
                    color=colors[char],
                    marker='o',
                    linestyle='',
                    markersize=5
                )
            ) for char in colors if char in detected_var
        ]
        ax.legend(
            custom_legend,
            [char for char in colors if char in detected_var],  # ensures correct sorting
            loc='lower right',
            title='variant',
            bbox_to_anchor=bbox_to_anchor,
            ncols=len(detected_var)/2 if aln.aln_type == 'AA' else len(detected_var),
            frameon=False
        )

    # format axis
    _format_x_axis(aln, ax, show_x_label, show_left=False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([ref_y_positions[x] for x in ref_y_positions])
    ax.set_yticklabels(ref_y_positions.keys())
    ax.set_ylim(0, y_pos)
    ax.set_ylabel('reference')


def orf_plot(aln: explore.MSA, ax: plt.Axes, min_length: int = 500, non_overlapping_orfs: bool = True, cmap: str = 'Blues', direction_marker_size: int | None = 5, show_x_label: bool = False, show_cbar: bool = False, cbar_fraction: float = 0.1):
    """
    Plot conserved ORFs.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param min_length: minimum length of orf
    :param non_overlapping_orfs: whether to consider overlapping orfs
    :param cmap: color mapping for % identity - see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    :param direction_marker_size: marker size for direction marker, not shown if marker_size == None
    :param show_x_label: whether to show the x-axis label
    :param show_cbar: whether to show the colorbar - see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
    :param cbar_fraction: fraction of the original ax reserved for the colorbar
    """

    # normalize colorbar
    cmap = ScalarMappable(norm=Normalize(0, 100), cmap=plt.get_cmap(cmap))

    # validate input
    _validate_input_parameters(aln, ax)

    # get orfs --> first deepcopy and reset zoom that the orfs are also zoomed in (by default, the orfs are only
    # searched within the zoomed region)
    aln_temp = deepcopy(aln)
    aln_temp.zoom = None
    if non_overlapping_orfs:
        annotation_dict = aln_temp.get_non_overlapping_conserved_orfs(min_length=min_length)
    else:
        annotation_dict = aln_temp.get_conserved_orfs(min_length=min_length)

    # filter dict for zoom
    if aln.zoom is not None:
        annotation_dict = {key:val for key, val in annotation_dict.items() if max(val['location'][0][0], aln.zoom[0]) <= min(val['location'][0][1], aln.zoom[1])}

    # add track for plotting
    _add_track_positions(annotation_dict)

    # plot
    _plot_annotation(annotation_dict, ax, direction_marker_size=direction_marker_size, color=cmap)

    # legend
    if show_cbar:
        cbar = plt.colorbar(cmap,ax=ax, location= 'top', orientation='horizontal', anchor=(1,0), shrink=0.2, pad=2/ax.bbox.height, fraction=cbar_fraction)
        cbar.set_label('% identity')
        cbar.set_ticks([0, 100])

    # format axis
    _format_x_axis(aln, ax, show_x_label, show_left=False)
    ax.set_ylim(bottom=0.8)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_title('conserved orfs', loc='left')


def annotation_plot(aln: explore.MSA, annotation: explore.Annotation | str, ax: plt.Axes, feature_to_plot: str, color: str = 'wheat', direction_marker_size: int | None = 5, show_x_label: bool = False):
    """
    Plot annotations from bed, gff or gb files. Are automatically mapped to alignment.
    :param aln: alignment MSA class
    :param annotation: annotation class | path to annotation file
    :param ax: matplotlib axes
    :param feature_to_plot: potential feature to plot (not for bed files as it is parsed as one feature)
    :param color: color for the annotation
    :param show_direction: show strand information
    :param direction_marker_size: marker size for direction marker, only relevant if show_direction is True
    :param show_x_label: whether to show the x-axis label
    """
    # helper function
    def parse_annotation_from_string(path: str, msa: explore.MSA) -> explore.Annotation:
        """
        Parse annotation.
        :param path: path to annotation
        :param msa: msa object
        :return: parsed annotation
        """
        if os.path.exists(path):
            # reset zoom so the annotation is correctly parsed
            msa_temp = deepcopy(msa)
            msa_temp.zoom = None
            return explore.Annotation(msa_temp, path)
        else:
            raise FileNotFoundError()

    # parse from path
    if type(annotation) is str:
        annotation = parse_annotation_from_string(annotation, aln)

    # validate input
    _validate_input_parameters(aln, ax, annotation)
    if not is_color_like(color):
        raise ValueError(f'{color} for reference is not a color')

    # ignore features to plot for bed files (here it is written into one feature)
    if annotation.ann_type == 'bed':
        annotation_dict = annotation.features['region']
        feature_to_plot = 'bed regions'
    else:
        # try to subset the annotation dict
        try:
            annotation_dict = annotation.features[feature_to_plot]
        except KeyError:
            raise KeyError(f'Feature {feature_to_plot} not found. Use annotation.features.keys() to see available features.')

    # plotting and formating
    _add_track_positions(annotation_dict)
    _plot_annotation(annotation_dict, ax, direction_marker_size=direction_marker_size, color=color)
    _format_x_axis(aln, ax, show_x_label, show_left=False)
    ax.set_ylim(bottom=0.8)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_title(f'{annotation.locus} ({feature_to_plot})', loc='left')

