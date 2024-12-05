

from msaexplorer import explore
from msaexplorer import config

from setuptools.errors import ClassError
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like
import numpy as np


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


def identity_plot(aln: explore.MSA, ax: plt.Axes, reference_color='blue', aln_colors:dict=config.ALN_COLORS, show_mask:bool=True, show_gaps:bool=True, show_mismatches:bool=True):
    """
    generates an alignment overview plot
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    :param reference_color: color of reference sequence
    :param aln_colors: dictionary containing colors: dict(0: color0, 1: color0, 2: color0, 3: color0) -> 0: match, 1: mismatch, 2: mask (N|X), 3: gap
    :param show_mask: whether to show N or X chars otherwise it will be shown as match or mismatch
    :param show_gaps: whether to show gaps otherwise it will be shown as match or mismatch
    :param show_mismatches: whether to show mismatches otherwise it will be shown as match
    """

    # input check
    if not isinstance(aln, explore.MSA):
        raise ClassError('alignment has to be an MSA class. use explore.MSA to read in alignment')
    if not isinstance(ax, plt.Axes):
        raise ClassError('ax has to be an matplotlib axis')

    # validate colors
    if not is_color_like(reference_color):
        raise ValueError('reference color is not a color')
    if aln_colors.keys() != config.ALN_COLORS.keys():
        raise ValueError('configure your dictionary with 0 to 3 key and associated colors. See config.ALN_COLORS')
    if not all([is_color_like(x) for x in aln_colors.values()]):
        raise ValueError('configure your dictionary with 0 to 3 key and associated colors. See config.ALN_COLORS')

    # get alignment and identity array
    alignment = aln.alignment
    identity_aln = aln.calc_identity_alignment(encode_mask=show_mask, encode_gaps=show_gaps, encode_mismatches=show_mismatches)
    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom
    # define the y position of the first sequence
    y_position = len(alignment) - 0.75
    # ini collection for patches
    col = []
    for sequence, seq_name in zip(identity_aln, alignment.keys()):
        # ini identity patches
        col.append(patches.Rectangle((zoom[0], y_position), zoom[1] - zoom[0],0.75,
                                               facecolor=aln_colors[0]
                                               )
                             )
        for identity_value in [3, 2, 1]: # first plot gaps, then mask, then mismatches
            stretches = find_stretches(sequence, identity_value)
            for stretch in stretches:
                col.append(
                    patches.Rectangle(
                        (stretch[0]+zoom[0], y_position),
                        stretch[1],
                        0.75,
                        color=aln_colors[identity_value],
                        linewidth=None
                    )
                )
        y_position -= 1

    ax.add_collection(PatchCollection(col, match_original=True, linewidths='none'))
    ax.set_ylim(0, len(alignment)+0.75/4)
    ax.set_xlim(zoom[0] - aln.length / 50, zoom[0] + aln.length + aln.length / 50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_xlabel('alignment position')