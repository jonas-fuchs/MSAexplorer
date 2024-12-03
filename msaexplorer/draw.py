from setuptools.errors import ClassError

from msaexplorer import explore
from msaexplorer import config

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np



def alignment_plot(aln: explore.MSA, ax: plt.Axes):
    """
    generates an alignment overview plot
    :param aln: alignment MSA class
    :param ax: matplotlib axes
    """

    # input check
    if not isinstance(aln, explore.MSA):
        raise ClassError('alignment has to be an MSA class. use explore.MSA to read in alignemnt')
    if not isinstance(ax, plt.Axes):
        raise ClassError('ax has to be an matplotlib axis')


    # get alignment and identity array
    alignment = aln.alignment
    identity_aln = aln.calc_identity_alignment()
    # define zoom to plot
    if aln.zoom is None:
        zoom = (0, aln.length)
    else:
        zoom = aln.zoom
    # define the y position of the first sequence
    y_position = len(alignment) - 0.75

    for sequence, seq_name in zip(identity_aln, alignment.keys()):

        # plot patch for the edge (excludes leading and trailing gaps)
        non_nan_indices = np.where(~np.isnan(sequence))[0]
        ax.add_patch(patches.Rectangle((int(non_nan_indices[0] + zoom[0]), y_position),
                                       int(non_nan_indices[-1] - non_nan_indices[0] + zoom[0]),
                                       0.75,
                                       facecolor='none',
                                       edgecolor='black',
                                       zorder=100)
                     )
        # ini indices
        previous_index = 0
        previous_position = sequence[0]
        for idx, position in enumerate(sequence):
            # for the last idx always plot
            if idx < aln.length - 1:
                # define plotting range of rectangles
                if position == previous_position or all(np.isnan([previous_position, position])):
                    continue
            if idx > previous_index:
                # define colors
                if np.isnan(previous_position):
                    facecolor = config.ALN_COLORS['del']
                elif alignment[seq_name][previous_index] in ['N', 'X']:
                    facecolor = config.ALN_COLORS['mask']
                else:
                    facecolor = config.ALN_COLORS[previous_position]
                # plot rectangle
                ax.add_patch(patches.Rectangle((previous_index + zoom[0], y_position),
                                               idx - previous_index + zoom[0],
                                               0.75,
                                               facecolor=facecolor)
                             )
                # additionally plot a vline --> ensures that also small differences in a large aln
                # are displayed
                ax.vlines(x=previous_index + zoom[0] + (idx - previous_index + zoom[0]) / 2,
                          ymin=y_position,
                          ymax=y_position + 0.75,
                          color=facecolor)
                # add a vline to gap sites (might be smaller in size as a single h line). likely needs a fix
                if facecolor == config.ALN_COLORS['del']:
                    ax.hlines(y=y_position + 0.75 / 2,
                              xmin=previous_index + zoom[0],
                              xmax=idx + zoom[0],
                              color='black')

            previous_index = idx
            previous_position = position

        y_position -= 1
    # adjust ax settings
    ax.set_ylim(0, len(alignment))
    ax.set_xlim(zoom[0] - aln.length / 50, zoom[0] + aln.length + aln.length / 50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.yticks([])
    plt.xlabel('alignment position')

