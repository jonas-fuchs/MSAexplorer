"""
contains the functions for drawing graphs
"""
# built-in
from typing import Callable, Dict
from copy import deepcopy

from numpy import ndarray

# MSAexplorer
from msaexplorer import explore
from msaexplorer import config

# libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.cm import ScalarMappable
from matplotlib.colors import is_color_like, Normalize, LinearSegmentedColormap
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# general helper functions
def _validate_input_parameters(aln: explore.MSA, ax: plt.Axes):
    """
    Validate MSA class and axis.
    """
    if not isinstance(aln, explore.MSA):
        raise ValueError('alignment has to be an MSA class. use explore.MSA() to read in alignment')
    if not isinstance(ax, plt.Axes):
        raise ValueError('ax has to be an matplotlib axis')


def _format_x_axis(aln: explore.MSA, ax: plt.Axes, show_x_label: bool, show_left: bool):
    """
    General axis formatting.
    """
    ax.set_xlim(
            (aln.zoom[0], aln.zoom[0] + aln.length) if aln.zoom is not None else (0, aln.length)
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
            ax.set_yticklabels(list(aln.alignment.keys())[::-1])
    else:
        ax.set_yticks([])


def _create_identity_patch(aln: explore.MSA, col: list, zoom: tuple[int, int], y_position: float | int, reference_color: str, seq_name: str, identity_color: str):
    """
    Creates the initial patch.
    """

    col.append(patches.Rectangle((zoom[0], y_position), zoom[1] - zoom[0],0.8,
                                                   facecolor=reference_color if seq_name == aln.reference_id else identity_color
                                                   )
                                 )


def _create_stretch_patch(col: list, stretches: list, zoom: tuple[int, int], y_position: float | int, colors: dict | ndarray, fancy_gaps: bool, matrix_value: int | float | ndarray):
    """
    Create a patch and add to list.
    """
    for stretch in stretches:
        col.append(
            patches.Rectangle(
                (stretch[0] + zoom[0], y_position),
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
                (stretch[0] + zoom[0], y_position + 0.375),
                stretch[1],
                0.05,
                color='black',
                linewidth=None
            )
        )


def identity_alignment(aln: explore.MSA, ax: plt.Axes, show_title: bool = True, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', aln_colors: dict = config.IDENTITY_COLORS, show_mask:bool = True, show_gaps:bool = True, fancy_gaps:bool = False, show_mismatches: bool = True, show_ambiguities: bool = False, show_x_label: bool = True, show_legend: bool = False, bbox_to_anchor: tuple[float|int, float|int] | list[float|int, float|int]= (1, 1.15)):
    """
    Generates an identity alignment overview plot.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param show_title: whether to show title
    :param show_seq_names: whether to show seq names
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param aln_colors: dictionary containing colors: dict(0: color0, 1: color0, 2: color0, 3: color0) -> 0: match, 1: mismatch, 2: mask (N|X), 3: gap
    :param show_mask: whether to show N or X chars otherwise it will be shown as match or mismatch
    :param show_gaps: whether to show gaps otherwise it will be shown as match or mismatch
    :param fancy_gaps: show gaps with a small black bar
    :param show_mismatches: whether to show mismatches otherwise it will be shown as match
    :param show_ambiguities: whether to show non-N ambiguities -> only relevant for RNA/DNA sequences
    :param show_x_label: whether to show x label
    :param show_legend: whether to show the legend
    :param bbox_to_anchor: bounding box coordinates for the legend - see: https://matplotlib.org/stable/api/legend_api.html
    """

    # input check
    _validate_input_parameters(aln, ax)

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError(f'{reference_color} for reference is not a color')
    if aln_colors.keys() != config.IDENTITY_COLORS.keys():
        raise ValueError('configure your dictionary like config.IDENTITY_COLORS')
    for key in aln_colors.keys():
        for key2 in aln_colors[key]:
            if key2 not in ['type', 'color']:
                raise ValueError('configure your dictionary like config.IDENTITY_COLORS')
    if not all([is_color_like(aln_colors[x]['color']) for x in aln_colors.keys()]):
        raise ValueError(f'one of the specified colors for the alignment is not a color')

    if fancy_gaps:
        show_gaps = True
    # get alignment and identity array
    identity_aln = aln.calc_identity_alignment(encode_mask=show_mask, encode_gaps=show_gaps, encode_mismatches=show_mismatches, encode_ambiguities=show_ambiguities)
    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom
    # define the y position of the first sequence
    y_position = len(aln.alignment) - 0.8
    detected_identity_values = {1}
    # ini collection for patches
    col = []
    for sequence, seq_name in zip(identity_aln, aln.alignment.keys()):
        # ini identity patches
        _create_identity_patch(aln, col, zoom, y_position, reference_color, seq_name, aln_colors[1]['color'])
        for identity_value in [np.nan, 0, 2, 3]:
            stretches = _find_stretches(sequence, identity_value)
            if not stretches:
                continue
            # add values for legend
            detected_identity_values.add(identity_value)
            _create_stretch_patch(col, stretches, zoom, y_position, aln_colors[identity_value]['color'], fancy_gaps, identity_value)
        y_position -= 1

    # custom legend
    if show_legend:
        custom_legend = [ax.add_line(plt.Line2D([], [], color=aln_colors[val]['color'], marker='s' ,markeredgecolor='grey', linestyle='', markersize=10)) for val in
                         detected_identity_values]
        ax.legend(
            custom_legend,
            [aln_colors[x]['type'] for x in detected_identity_values],
            loc='upper right',
            bbox_to_anchor=bbox_to_anchor,
            ncols=len(detected_identity_values),
            frameon=False
        )
    # format seq names
    _seq_names(aln, ax, custom_seq_names, show_seq_names)
    # configure axis
    ax.add_collection(PatchCollection(col, match_original=True, linewidths='none', joinstyle='miter', capstyle='butt'))
    ax.set_ylim(0, len(aln.alignment)+0.2)
    if show_title:
        ax.set_title('identity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def similarity_alignment(aln: explore.MSA, ax: plt.Axes, matrix_type: str | None = None, show_title: bool = True, show_seq_names: bool = False, custom_seq_names: tuple | list = (), reference_color: str = 'lightsteelblue', similarity_colors: tuple[str, str] | list[str, str] = ('darkblue', 'lightgrey'), gap_color: str = 'white', show_gaps:bool = True, fancy_gaps:bool = False, show_x_label: bool = True, show_cbar: bool = False, cbar_fraction: float = 0.1):
    """
    Generates a similarity alignment overview plot.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param matrix_type: substitution matrix - see config.SUBS_MATRICES, standard: NT - TRANS, AA - BLOSUM65
    :param show_title: whether to show title
    :param show_seq_names: whether to show seq names
    :param custom_seq_names: custom seq names
    :param reference_color: color of reference sequence
    :param similarity_colors: first color - not similar, second color - similar --> create a colormap
    :param gap_color: color for gaps
    :param show_gaps: whether to show gaps otherwise it will be ignored
    :param fancy_gaps: show gaps with a small black bar
    :param show_x_label: whether to show x label
    :param show_cbar: whether to show the legend - see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
    :param cbar_fraction: fraction of the original ax reserved for the legend
    """
    # input check
    _validate_input_parameters(aln, ax)
    # always show gaps for fancy gaps
    if fancy_gaps:
        show_gaps = True
    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError(f'{reference_color} for reference is not a color')
    for color in similarity_colors:
        if not is_color_like(color):
            raise ValueError(f'{color} for similarity is not a color')

    # get data
    similarity_aln = aln.calc_similarity_alignment(matrix_type=matrix_type)

    # determine min max values of the underlying matrix and create cmap
    matrix = config.SUBS_MATRICES[aln.aln_type][similarity_aln.dtype.metadata['matrix']]
    min_value, max_value = min([min(matrix[d].values()) for d in matrix]), max([max(matrix[d].values()) for d in matrix])
    cmap = ScalarMappable(norm=Normalize(
        vmin=min_value,
        vmax=max_value
    ), cmap=LinearSegmentedColormap.from_list('', colors=similarity_colors))
    # create plot
    col = []
    y_position = len(aln.alignment) - 0.8
    for sequence, seq_name in zip(similarity_aln, aln.alignment.keys()):
        _create_identity_patch(aln, col, zoom, y_position, reference_color, seq_name, similarity_colors[1])
        similarity_values = np.append([np.nan], np.arange(start=min_value, stop=max_value), 0) if show_gaps else np.arange(start=min_value, stop=max_value)
        for similarity_value in similarity_values:
            stretches = _find_stretches(sequence, similarity_value)
            if not stretches:
                continue
            _create_stretch_patch(col, stretches, zoom, y_position, cmap.to_rgba(similarity_value) if not np.isnan(similarity_value) else gap_color,
                                  fancy_gaps, similarity_value)
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
    ax.set_ylim(0, len(aln.alignment)+0.2)
    if show_title:
        ax.set_title('similarity', loc='left')
    _format_x_axis(aln, ax, show_x_label, show_left=False)


def stat_plot(aln: explore.MSA, ax: plt.Axes, stat_type: str, line_color: str = 'burlywood', line_width: int | float = 2, rolling_average: int = 20, show_x_label: bool = False):
    """
    Generate a plot for the various alignment stats.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param stat_type: 'entropy', 'gc', 'coverage', 'identity' or 'similarity' -> (here default matrices are used NT - TRANS, AA - BLOSUM65)
    :param line_color: color of the line
    :param line_width: width of the line
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
    stat_functions: Dict[str, Callable[[], list | ndarray]] = {
        'gc': aln.calc_gc,
        'entropy': aln.calc_entropy,
        'coverage': aln.calc_coverage,
        'identity': aln.calc_identity_alignment,
        'similarity': aln.calc_similarity_alignment,
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
    array = stat_functions[stat_type]()
    if stat_type == 'similarity':
        matrix = config.SUBS_MATRICES[aln.aln_type][array.dtype.metadata['matrix']]
        min_value, max_value = min([min(matrix[d].values()) for d in matrix]), max([max(matrix[d].values()) for d in matrix])
    else:
        min_value, max_value = 0, 1

    if stat_type in ['identity', 'similarity']:
        # for the mean nan values get handled as the lowest possible number in the matrix
        array = np.nan_to_num(array, True, 0 if stat_type == 'identity' else min_value)
        array = np.mean(array, axis=0)

    data, plot_idx = moving_average(array, rolling_average)
    # plot
    ax.plot(plot_idx, data, color=line_color, linewidth=line_width)

    # specific visual cues for individual plots
    if stat_type != 'gc':
        ax.fill_between(plot_idx, data, color=(line_color, 0.5))
    else:
        ax.hlines(0.5, xmin=0, xmax=aln.zoom[0] + aln.length if aln.zoom is not None else aln.length, color='black', linestyles='--', linewidth=1)

    # format axis
    ax.set_ylim(min_value, max_value*0.1+max_value)
    ax.set_yticks([min_value, max_value])
    if stat_type != 'gc':
        ax.set_yticklabels(['low', 'high'])
    else:
        ax.set_yticklabels(['0', '100'])

    _format_x_axis(aln, ax, show_x_label, show_left=True)
    if rolling_average > 1:
        ax.set_ylabel(f'{stat_type}\n {2*rolling_average}c mean')
    else:
        ax.set_ylabel(f'{stat_type}')


def variant_plot(aln: explore.MSA, ax: plt.Axes, lollisize: tuple[int, int] | list[int, int] = (1, 3), show_x_label: bool = False, colors: dict | None = None, show_legend: bool = True, bbox_to_anchor: tuple[float|int, float|int] | list[float|int, float|int] = (1, 1.15)):
    """
    Plots variants.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param lollisize: (stem_size, head_size)
    :param show_x_label:  whether to show the x-axis label
    :param colors: colors for variants - if None standard colors are used (config.AA_colors or config.NT_colors)
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
        for char in colors:
            if not is_color_like(colors[char]):
                raise TypeError(f'{colors[char]} is not a color')

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
            bbox_to_anchor=bbox_to_anchor,
            ncols=len(detected_var),
            frameon=False
        )

    # format axis
    _format_x_axis(aln, ax, show_x_label, show_left=False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([ref_y_positions[x] for x in ref_y_positions])
    ax.set_yticklabels(ref_y_positions.keys())
    ax.set_ylim(0, y_pos)
    ax.set_ylabel('reference')


def orf_plot(aln: explore.MSA, ax: plt.Axes, min_length: int = 500, non_overlapping_orfs: bool = True, cmap: str = 'Blues', show_direction:bool = True, direction_marker_size: int = 5, show_x_label: bool = False, show_cbar: bool = False, cbar_fraction: float = 0.1):
    """
    Plot conserved ORFs.
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param min_length: minimum length of orf
    :param non_overlapping_orfs: whether to consider overlapping orfs
    :param cmap: color mapping for % identity - see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    :param show_direction: show strand information
    :param direction_marker_size: marker size for direction marker, only relevant if show_direction is True
    :param show_x_label: whether to show the x-axis label
    :param show_cbar: whether to show the colorbar - see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
    :param cbar_fraction: fraction of the original ax reserved for the colorbar

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
        annotation_dict = {key:val for key, val in annotation_dict.items() if aln.zoom[0] < val['positions'][0] <= aln.zoom[1]}

    add_track_positions(annotation_dict)

    # plot
    max_track = 0
    for annotation in annotation_dict:
        x_value = annotation_dict[annotation]['positions'][0]
        length = annotation_dict[annotation]['positions'][1] - annotation_dict[annotation]['positions'][0]
        ax.add_patch(
            patches.FancyBboxPatch(
                (x_value, annotation_dict[annotation]['track'] + 1),
                length,
                0.8,
                boxstyle="round",
                ec="black",
                fc=cmap.to_rgba(annotation_dict[annotation]['conservation'])
            )
        )
        if annotation_dict[annotation]['track'] > max_track:
            max_track = annotation_dict[annotation]['track']

        if show_direction:
            if annotation_dict[annotation]['strand'] == '-':
                marker = '<'
            else:
                marker = '>'
            ax.plot(x_value + length/2, annotation_dict[annotation]['track'] + 1.4, marker=marker, markersize=direction_marker_size, color='white', markeredgecolor='black')

    # legend
    if show_cbar:
        cbar = plt.colorbar(cmap,ax=ax, location= 'top', orientation='horizontal', anchor=(1,0), shrink=0.2, pad=2/ax.bbox.height, fraction=cbar_fraction)
        cbar.set_label('% identity')
        cbar.set_ticks([0, 100])

    # format axis
    _format_x_axis(aln, ax, show_x_label, show_left=False)
    ax.set_ylim(bottom=0)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_title('conserved orfs', loc='left')

