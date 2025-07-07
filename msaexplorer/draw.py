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
from matplotlib.collections import PatchCollection, PolyCollection


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


def _find_stretches(row, non_nan_only=False) -> list[tuple[int, int, int]] | list[tuple[int, int]] :
    """
    Finds consecutive stretches of values in an array, with an option to exclude NaN stretches.

    :param row: NumPy array of values
    :param non_nan_only: If True, only returns stretches of non-NaN values.
    :return: List of stretches (start, end, value at start) for all values or non-NaN values only.
    """
    if row.size == 0:
        return []

    if non_nan_only:
        # Create a boolean mask for non-NaN values
        non_nan_mask = ~np.isnan(row)
        # Find changes in the mask
        changes = np.diff(non_nan_mask.astype(int)) != 0
        change_idx = np.nonzero(changes)[0]
        starts = np.concatenate(([0], change_idx + 1))
        ends = np.concatenate((change_idx, [len(row) - 1]))

        # Return only stretches that start with non-NaN values
        return [(start, end) for start, end in zip(starts, ends) if non_nan_mask[start]]
    else:
        # Find change points: where adjacent cells differ.
        changes = np.diff(row) != 0
        change_idx = np.nonzero(changes)[0]
        starts = np.concatenate(([0], change_idx + 1))
        ends = np.concatenate((change_idx, [len(row) - 1]))

        return [(start, end, row[start]) for start, end in zip(starts, ends)]


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
        ax.set_yticks(np.arange(len(aln.alignment)) + 0.4)
        if custom_seq_names:
            ax.set_yticklabels(custom_seq_names[::-1])
        else:
            names = [x.split(' ')[0] for x in list(aln.alignment.keys())[::-1]]
            ax.set_yticklabels(names)
    else:
        ax.set_yticks([])


def _create_identity_patch(row, aln: explore.MSA, col: list, zoom: tuple[int, int], y_position: float | int, reference_color: str, seq_name: str, identity_color: str | ndarray, show_gaps: bool):
    """
    Creates the initial patch.
    """
    if show_gaps:
        # plot a rectangle for parts that do not have gaps
        for stretch in _find_stretches(row, True):
            col.append(patches.Rectangle((stretch[0] + zoom[0] - 0.5, y_position), stretch[1] - stretch[0] + 1, 0.8,
                                                           facecolor=reference_color if seq_name == aln.reference_id else identity_color
                                                           )

                                         )
    # just plot a rectangle
    else:
        col.append(patches.Rectangle((zoom[0] - 0.5, y_position), zoom[1] - zoom[0], 0.8,
                                     facecolor=reference_color if seq_name == aln.reference_id else identity_color
                                     )
                   )


def _create_polygons(stretches: list, identity_values: list | ndarray, zoom: tuple, y_position: int, polygons: list, aln_colors: dict | ScalarMappable, polygon_colors: list, detected_identity_values: set | None = None):
    """
    create the individual polygons
    """

    for start, end, value in stretches:
        if value not in identity_values:
            continue
        if detected_identity_values is not None:
            detected_identity_values.add(value)
        width = end + 1 - start
        # Calculate x coordinates adjusted for zoom and centering
        x0 = start + zoom[0] - 0.5
        x1 = x0 + width
        # Define the rectangle corners
        rect_coords = [
            (x0, y_position),
            (x1, y_position),
            (x1, y_position + 0.8),
            (x0, y_position + 0.8),
            (x0, y_position)
        ]
        polygons.append(rect_coords)
        if type(aln_colors) != ScalarMappable:
            polygon_colors.append(aln_colors[value]['color'])
        else:
            polygon_colors.append(aln_colors.to_rgba(value))


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
    return 'white' if brightness < 0.5 else 'black'


def _plot_sequence_text(aln: explore.MSA, seq_name: str, ref_name: str, allways_text: bool, values: ndarray, matrix: ndarray, ax:plt.Axes, zoom: tuple, y_position:int, value_to_skip: int, ref_color:str, show_gaps: bool, cmap: None | ScalarMappable = None):
    """
    Plot sequence text - however this will be done even if there is not enough space.
    Might need some rework in the future.
    """
    x_text = 0
    if seq_name == ref_name:
        different_cols = np.any((matrix != value_to_skip) & ~np.isnan(matrix), axis=0)
    else:
        different_cols = [False]*aln.length

    for idx, (character, value) in enumerate(zip(aln.alignment[seq_name], values)):
        if value != value_to_skip and character != '-' or seq_name == ref_name and character != '-' or character == '-'and not show_gaps or allways_text and character != '-':
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
                c=text_color,
            )
        x_text += 1


def identity_alignment(aln: explore.MSA, ax: plt.Axes, show_title: bool = True, show_identity_sequence: bool = False, show_sequence_all: bool = False, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', show_mask:bool = True, show_gaps:bool = True, fancy_gaps:bool = False, show_mismatches: bool = True, show_ambiguities: bool = False, color_scheme: str | None = None, show_x_label: bool = True, show_legend: bool = False, bbox_to_anchor: tuple[float|int, float|int] | list[float|int, float|int]= (1, 1)):
    """
    Generates an identity alignment overview plot.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param show_title: whether to show title
    :param show_seq_names: whether to show seq names
    :param show_identity_sequence: whether to show sequence for only differences and reference - zoom in to avoid plotting issues
    :param show_sequence_all: whether to show all sequences - zoom in to avoid plotting issues
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param show_mask: whether to show N or X chars otherwise it will be shown as match or mismatch
    :param show_gaps: whether to show gaps otherwise it will be shown as match or mismatch
    :param fancy_gaps: show gaps with a small black bar
    :param show_mismatches: whether to show mismatches otherwise it will be shown as match
    :param show_ambiguities: whether to show non-N ambiguities -> only relevant for RNA/DNA sequences
    :param color_scheme: color mismatching chars with their unique color. Options for DNA/RNA are: standard, purine_pyrimidine, strong_weak; and for AS: standard, clustal, zappo, hydrophobicity
    :param show_x_label: whether to show x label
    :param show_legend: whether to show the legend
    :param bbox_to_anchor: bounding box coordinates for the legend - see: https://matplotlib.org/stable/api/legend_api.html
    """

    # Validate inputs and colors
    _validate_input_parameters(aln, ax)
    if not is_color_like(reference_color):
        raise ValueError(f'{reference_color} for reference is not a color')

    if color_scheme is not None:
        if aln.aln_type == 'AA' and color_scheme not in ['standard', 'clustal', 'zappo', 'hydrophobicity']:
            raise ValueError(f'{color_scheme} is not a supported coloring scheme for {aln.aln_type} alignments. Supported are: standard, clustal, zappo, hydrophobicity')
        elif aln.aln_type in ['DNA', 'RNA'] and color_scheme not in ['standard', 'purine_pyrimidine', 'strong_weak']:
            raise ValueError(f'{color_scheme} is not a supported coloring scheme for {aln.aln_type} alignments. Supported are: standard, purine_pyrimidine, strong_weak')

    # Both options for gaps work hand in hand
    if fancy_gaps:
        show_gaps = True

    # Determine zoom
    zoom = (0, aln.length) if aln.zoom is None else aln.zoom

    # Set up color mapping and identity values
    aln_colors = config.IDENTITY_COLORS.copy()
    identity_values = [-1, -2, -3]  # -1 = mismatch, -2 = mask, -3 ambiguity
    if color_scheme is not None:
        colors_to_extend = config.CHAR_COLORS[aln.aln_type][color_scheme]
        identity_values += [x + 1 for x in range(len(colors_to_extend))] # x+1 is needed to allow correct mapping
        # use the standard setting for the index (same as in aln.calc_identity_alignment)
        # and map the corresponding color scheme to it
        for idx, char in enumerate(config.CHAR_COLORS[aln.aln_type]['standard']):
            aln_colors[idx + 1] = {'type': char, 'color': colors_to_extend[char]}

    # Compute identity alignment
    identity_aln = aln.calc_identity_alignment(
        encode_mask=show_mask,
        encode_gaps=show_gaps,
        encode_mismatches=show_mismatches,
        encode_ambiguities=show_ambiguities,
        encode_each_mismatch_char=True if color_scheme is not None else False
    )

    # List to store polygons
    detected_identity_values = {0}
    polygons, polygon_colors, patch_list = [], [], []

    for i, seq_name in enumerate(aln.alignment):
        y_position = len(aln.alignment) - i - 1
        row = identity_aln[i]

        # plot a line below everything
        if fancy_gaps:
            ax.hlines(
                y_position + 0.4,
                xmin=zoom[0] - 0.5,
                xmax=zoom[1] + 0.5,
                color='black',
                linestyle='-',
                zorder=0,
                linewidth=0.75
            )

        # plot the basic shape per sequence with gaps
        _create_identity_patch(row, aln, patch_list, zoom, y_position, reference_color, seq_name, aln_colors[0]['color'], show_gaps)

        # find consecutive stretches
        stretches = _find_stretches(row)
        # create polygons per stretch
        _create_polygons(stretches, identity_values, zoom, y_position, polygons, aln_colors, polygon_colors, detected_identity_values)

        # add sequence text
        if show_identity_sequence or show_sequence_all:
            _plot_sequence_text(aln, list(aln.alignment.keys())[i], aln.reference_id, show_sequence_all, identity_aln[i], identity_aln, ax,
                                zoom, y_position, 0, reference_color, show_gaps)

    # Create the LineCollection: each segment is drawn in a single call.
    ax.add_collection(PatchCollection(patch_list, match_original=True, linewidths='none', joinstyle='miter', capstyle='butt'))
    ax.add_collection(PolyCollection(polygons, facecolors=polygon_colors, linewidths=0.5, edgecolors=polygon_colors))

    # custom legend
    if show_legend:
        if color_scheme is not None and color_scheme != 'standard':
            for x in aln_colors:
                for group in config.CHAR_GROUPS[aln.aln_type][color_scheme]:
                    if aln_colors[x]['type'] in config.CHAR_GROUPS[aln.aln_type][color_scheme][group]:
                        aln_colors[x]['type'] = group
                        break
        # create it
        handels, labels, detected_groups = [], [], set()
        for x in aln_colors:
            if x in detected_identity_values and aln_colors[x]['type'] not in detected_groups:
                handels.append(
                    ax.add_line(
                        plt.Line2D(
                            [],
                            [],
                            color=aln_colors[x]['color'], marker='s', markeredgecolor='grey', linestyle='',
                            markersize=10))
                )
                labels.append(aln_colors[x]['type'])
                detected_groups.add(aln_colors[x]['type'])

        # plot it
        ax.legend(
            handels,
            labels,
            loc='lower right',
            bbox_to_anchor=bbox_to_anchor,
            ncols=(len(detected_identity_values) + 1) / 2 if aln.aln_type == 'AA' and color_scheme == 'standard' else len(detected_identity_values),
            frameon=False
        )

    _seq_names(aln, ax, custom_seq_names, show_seq_names)

    # configure axis
    ax.set_ylim(0, len(aln.alignment))
    if show_title:
        ax.set_title('identity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def similarity_alignment(aln: explore.MSA, ax: plt.Axes, matrix_type: str | None = None, show_title: bool = True, show_identity_sequence: bool = False, show_sequence_all: bool = False, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', cmap: str = 'twilight_r', show_gaps:bool = True, fancy_gaps:bool = False, show_x_label: bool = True, show_cbar: bool = False, cbar_fraction: float = 0.1):
    """
    Generates a similarity alignment overview plot. Importantly the similarity values are normalized!
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param matrix_type: substitution matrix - see config.SUBS_MATRICES, standard: NT - TRANS, AA - BLOSUM65
    :param show_title: whether to show title
    :param show_identity_sequence: whether to show sequence only for differences and reference - zoom in to avoid plotting issues
    :param show_sequence_all: whether to show all sequences - zoom in to avoid plotting issues
    :param show_seq_names: whether to show seq names
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param cmap: color mapping for % identity - see https://matplotlib.org/stable/users/explain/colors/colormaps.html
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
    similarity_values = np.arange(start=min_value, stop=max_value, step=0.01)
    # round it to be absolutely sure that values match with rounded sim alignment
    similarity_values = similarity_values.round(2)

    # create plot
    polygons, polygon_colors, patch_list = [], [], []

    for i, seq_name in enumerate(aln.alignment):
        y_position = len(aln.alignment) - i - 1
        row = similarity_aln[i]

        # plot a line below everything
        if fancy_gaps:
            ax.hlines(
                y_position + 0.4,
                xmin=zoom[0] - 0.5,
                xmax=zoom[1] + 0.5,
                color='black',
                linestyle='-',
                zorder=0,
                linewidth=0.75
            )

        # plot the basic shape per sequence with gaps
        _create_identity_patch(row, aln, patch_list, zoom, y_position, reference_color, seq_name,
                               cmap.to_rgba(max_value) if seq_name != aln.reference_id else reference_color,
                               show_gaps)

        # find consecutive stretches
        stretches = _find_stretches(row)
        # create polygons per stretch
        _create_polygons(stretches, similarity_values, zoom, y_position, polygons, cmap, polygon_colors)

        # add sequence text
        if show_sequence_all or show_identity_sequence:
            _plot_sequence_text(aln, list(aln.alignment.keys())[i], aln.reference_id, show_sequence_all, similarity_aln[i], similarity_aln, ax,
                                zoom, y_position, 1, reference_color, show_gaps, cmap)

    # Create the LineCollection: each segment is drawn in a single call.
    ax.add_collection(PatchCollection(patch_list, match_original=True, linewidths='none', joinstyle='miter', capstyle='butt'))
    ax.add_collection(PolyCollection(polygons, facecolors=polygon_colors, linewidths=0.5, edgecolors=polygon_colors))

    # legend
    if show_cbar:
        cbar = plt.colorbar(cmap, ax=ax, location= 'top', anchor=(1,0), shrink=0.2, pad=2/ax.bbox.height, fraction=cbar_fraction)
        cbar.set_ticks([min_value, max_value])
        cbar.set_ticklabels(['low', 'high'])

    # format seq names
    _seq_names(aln, ax, custom_seq_names, show_seq_names)

    # configure axis
    ax.set_ylim(0, len(aln.alignment))
    if show_title:
        ax.set_title('similarity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def _moving_average(arr: ndarray, window_size: int, zoom: tuple | None, aln_length: int) -> tuple[ndarray, ndarray]:
    """
    Calculate the moving average of an array.
    :param arr: array with values
    :param window_size: size of the moving average
    :param zoom: zoom of the alignment
    :param aln_length: length of the alignment
    :return: new array with moving average
    """
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

        return np.array(moving_averages), np.array(plotting_idx) if zoom is None else np.array(plotting_idx) + zoom[0]
    else:
        return arr, np.arange(zoom[0], zoom[1]) if zoom is not None else np.arange(aln_length)


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
    data, plot_idx = _moving_average(array, rolling_average, aln.zoom, aln.length)

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
    colors = config.CHAR_COLORS[aln.aln_type]['standard']
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

#TODO plot a stacked bar/sequence logo based on PWM