"""
This contains the main functions to create the plots in the app
"""

# libs
import math
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from msaexplorer import draw

def set_aln(aln, inputs):
    # set the reference sequence
    if 'reference' in inputs:
        if inputs['reference'] == 'first':
            aln.reference_id = list(aln.alignment.keys())[0]
        elif inputs['reference'] == 'consensus':
            aln.reference_id = None
        else:
            aln.reference_id = inputs['reference']

    # Update zoom level from slider -> +1 needed as msaexplorer uses range for zoom
    if 'zoom_range' in inputs:
        aln.zoom = (inputs['zoom_range'][0], inputs['zoom_range'][1] + 1)

    return aln


def create_msa_plot(aln, ann, inputs, fig_size=None) -> plt.Figure | None:
    """
    :param aln: MSA object
    :param ann: Annotation object
    :param inputs: all user inputs
    :param fig_size: size of the figure -> for pdf plotting
    :return: figure
    """
    if not aln:
        return None

    aln = set_aln(aln, inputs)

    # Collect height ratios and corresponding axes
    height_ratios = []
    plot_functions = []

    # First plot
    if inputs['stat_type'] not in ['Off', 'sequence logo']:
        height_ratios.append(inputs['plot_1_size'])
        plot_functions.append(
            lambda ax: draw.stat_plot(
                aln,
                ax=ax,
                stat_type=inputs['stat_type'],
                line_width=1,
                rolling_average=inputs['rolling_average'],
                show_x_label=True if inputs['annotation'] == 'Off' and inputs['alignment_type'] == 'Off' else False,
                line_color=inputs['stat_color'])
        )
    elif inputs['stat_type'] == 'sequence logo':
        height_ratios.append(inputs['plot_1_size'])
        plot_functions.append(
            lambda ax: draw.sequence_logo(
                aln, ax=ax,
                show_x_label=True if inputs['annotation'] == 'Off' and inputs['alignment_type'] == 'Off' else False,
                plot_type = inputs['logo_type'],
                color_scheme = inputs['logo_coloring']
            )
        )

    # Second plot
    if inputs['alignment_type'] != 'Off':
        height_ratios.append(inputs['plot_2_size'])
        plot_functions.append(
            lambda ax: draw.identity_alignment(
                aln, ax=ax,
                show_identity_sequence=inputs['show_sequence'],
                show_sequence_all=inputs['show_sequence_all'],
                fancy_gaps=inputs['fancy_gaps'],
                show_gaps=inputs['show_gaps'],
                show_mask=inputs['show_mask'],
                show_mismatches=True,
                show_ambiguities=inputs['show_ambiguities'],
                color_scheme=inputs['char_coloring'] if inputs['char_coloring'] != 'None' else None,
                identical_char_color=inputs['identical_char_color'],
                different_char_color=inputs['different_char_color'],
                mask_color=inputs['mask_color'],
                ambiguity_color=inputs['ambiguity_color'],
                reference_color=inputs['reference_color'],
                show_seq_names=inputs['seq_names'],
                show_x_label=True if inputs['annotation'] == 'Off' else False,
                show_legend=inputs['show_legend']
        ) if inputs['alignment_type'] == 'identity' else draw.alignment(
                aln, ax=ax,
                show_sequence_all=inputs['show_sequence_all'],
                fancy_gaps=inputs['fancy_gaps'],
                show_mask=inputs['show_mask'],
                show_ambiguities=inputs['show_ambiguities'],
                color_scheme=inputs['char_coloring'] if inputs['char_coloring'] != 'None' else None,
                mask_color=inputs['mask_color'],
                ambiguity_color=inputs['ambiguity_color'],
                show_seq_names=inputs['seq_names'],
                show_x_label=True if inputs['annotation'] == 'Off' else False,
                show_legend=inputs['show_legend']
        ) if inputs['alignment_type'] == 'normal' else draw.similarity_alignment(
                aln, ax=ax,
                show_similarity_sequence=inputs['show_sequence'],
                show_sequence_all=inputs['show_sequence_all'],
                fancy_gaps=inputs['fancy_gaps'],
                show_gaps=inputs['show_gaps'],
                reference_color=inputs['reference_color'],
                matrix_type=inputs['matrix'],
                show_seq_names=inputs['seq_names'],
                cmap=inputs['matrix_color_mapping'],
                show_cbar=inputs['show_legend'],
                cbar_fraction=0.02,
                show_x_label=True if inputs['annotation'] == 'Off' else False)
    )

    # Third Plot
    if inputs['annotation'] != 'Off':
        height_ratios.append(inputs['plot_3_size'])
        plot_functions.append(
            lambda ax: draw.annotation_plot(
                aln, ann, ax=ax,
                feature_to_plot=inputs['feature_display'],
                color=inputs['feature_color'],
                direction_marker_size=inputs['strand_marker_size'],
                show_x_label=True
        ) if (inputs['annotation'] == 'Annotation' and inputs['annotation_file']) or (inputs['annotation'] == 'Annotation' and inputs['example_alignment']) else draw.orf_plot(
                aln, ax=ax,
                cmap=inputs['color_mapping'],
                non_overlapping_orfs=inputs['non_overlapping'],
                direction_marker_size=inputs['strand_marker_size'],
                show_x_label=True,
                show_cbar=inputs['show_legend_third_plot'],
                cbar_fraction=0.2,
                min_length=inputs['min_orf_length']
        ) if inputs['annotation'] == 'Conserved ORFs' else draw.variant_plot(
                aln, ax=ax,
                show_x_label=True,
                lollisize=(inputs['stem_size'], inputs['head_size']),
                color_scheme=inputs['snp_coloring'],
                show_legend=inputs['show_legend_third_plot'])
    )
    # do not plot anything if all plots are off
    if not height_ratios:
        return None

    # Prepare the plot with dynamic number of subplots
    if fig_size is not None:
        fig, axes = plt.subplots(nrows=len(height_ratios), height_ratios=height_ratios, figsize=fig_size)
    else:
        fig, axes = plt.subplots(nrows=len(height_ratios), height_ratios=height_ratios)

    # If there is only one plot, `axes` is not a list
    if len(height_ratios) == 1:
        axes = [axes]

    # Render each enabled plot
    for ax, plot_func in zip(axes, plot_functions):
        plot_func(ax)

    return fig


def create_analysis_custom_heatmap(aln, inputs):
    """
    Create a heatmap that depends on the user (multiple options)
    """
    if aln is None:
        return None

    # scale the plot for 70 % of the window
    if inputs['dimensions']['width'] > inputs['dimensions']['height']:
        figure_size = int(inputs['dimensions']['height'] * 0.7)
    else:
        figure_size = int(inputs['dimensions']['width'] * 0.7)

    matrix = aln.calc_pairwise_identity_matrix(inputs['additional_analysis_options_left'])
    labels = [x.split(' ')[0] for x in list(aln.alignment.keys())]

    # generate hover text
    hover_text = [
        [
            f'{labels[i]} & {labels[j]}: {matrix[i][j]:.2f}%'
            for j in range(len(matrix[i]))
        ]
        for i in range(len(matrix))
    ]

   # Create the heatmap
    fig = go.Figure(
        go.Heatmap(
            z=matrix,
            x=labels,
            y=labels,
            text=hover_text,
            hoverinfo='text',
            colorscale='deep',
            colorbar=dict(thickness=20, ticklen=4, len=0.8, title='%')
        )
    )

    fig.update_layout(
        title='Pairwise identities',
        width=figure_size,
        height=figure_size - 50,
        margin=dict(t=50, b=50, l=50, r=50),
        xaxis=dict(scaleanchor='y', constrain='domain'),
        yaxis=dict(scaleanchor='x', constrain='domain')
    )

    return fig


def create_freq_heatmap(aln, inputs):
    """
    Create a frequency heatmap
    """
    if aln is None:
        return None

    # analog to the other heatmap
    if inputs['dimensions']['width'] > inputs['dimensions']['height']:
        figure_height = int(inputs['dimensions']['height'] * 0.7)
    else:
        figure_height = int(inputs['dimensions']['width'] * 0.7)

    # get char freqs
    freqs = aln.calc_character_frequencies()
    characters = sorted(set(char for seq_id in freqs if seq_id != 'total' for char in freqs[seq_id] if char not in ['-', 'N', 'X']))
    sequence_ids = [seq_id for seq_id in freqs if seq_id != 'total']

    # create the matrix
    matrix, seq_ids = [], []
    for seq_id in sequence_ids:
        row = []
        seq_ids.append(seq_id.split(' ')[0])
        for char in characters:
            row.append(freqs[seq_id].get(char, {'% of non-gapped': 0})['% of non-gapped'])
        matrix.append(row)

    matrix = np.array(matrix)

    # generate hover text
    hover_text = [
        [
            f'{seq_ids[i]}: {matrix[i][j]:.2f}% {characters[j]}'
            for j in range(len(matrix[i]))
        ]
        for i in range(len(matrix))
    ]

    # plot
    fig = go.Figure(
            go.Heatmap(
                z=matrix,
                x=characters,
                y=seq_ids,
                text=hover_text,
                hoverinfo='text',
                colorscale='Cividis',
                colorbar=dict(thickness=20, ticklen=4, title='%', len=0.7),
            )
        )

    fig.update_layout(
        title='Character frequencies (% ungapped)',
        height=figure_height - 50,
        margin=dict(t=50, b=50, l=50, r=50),
    )

    return fig


def create_recovery_heatmap(aln, inputs):
    """
    Create a frequency heatmap (right heatmap)
    """

    if aln is None:
        return None

    # analog to the other heatmap
    if inputs['dimensions']['width'] > inputs['dimensions']['height']:
        figure_height = int(inputs['dimensions']['height'] * 0.7)
    else:
        figure_height = int(inputs['dimensions']['width'] * 0.7)

    # ini
    recovery_dict = aln.calc_percent_recovery()
    sequence_ids = [x.split()[0] for x in recovery_dict.keys()]
    recovery_vals = list(recovery_dict.values())
    n = len(sequence_ids)
    side = math.ceil(math.sqrt(n))

    # fill the grid
    padded_length = side * side
    sequence_ids += [0] * (padded_length - n)
    recovery_vals += [0] * (padded_length - n)

    # get coordinates
    xs = [i % side for i in range(padded_length)]
    ys = [i // side for i in range(padded_length)]

    # plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs,
        y=ys,
        mode='markers',
        marker=dict(
            size=recovery_vals,
            sizemode='diameter',
            sizeref=2.5,
            sizemin=4,
            cmin=0,
            cmax=100,
            color=recovery_vals,
            colorscale='RdBu',
            colorbar=dict(thickness=20, ticklen=4, title='%', len=0.7),
            line=dict(width=0)
        ),
        text=[
            f"{sid}<br>{round(val, 0)}% recovery" if sid is not None and val !=0 else None
            for sid, val in zip(sequence_ids, recovery_vals)
        ],
        hoverinfo='text'
    ))

    ref = aln.reference_id.split()[0] if aln.reference_id is not None else 'consensus'

    fig.update_layout(
        title=f"Recovery compared to: {ref}",
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        plot_bgcolor="white",
        height=figure_height - 50,
        margin=dict(t=50, b=50, l=50, r=50),
    )

    return fig