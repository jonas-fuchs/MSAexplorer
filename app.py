"""
This contains the code to create the MSAexplorer shiny application
"""
# build-in
from pathlib import Path
import tempfile

import numpy as np
# libs
from shiny import App, render, ui, reactive
from shinywidgets import output_widget, render_widget
import plotly.graph_objects as go
import matplotlib
from matplotlib import colormaps
import matplotlib.pyplot as plt

# msaexplorer
from msaexplorer import explore, draw, config, export

# file paths for css and js
css_file = Path(__file__).resolve().parent / 'www' /'css' / 'styles.css'
js_file = Path(__file__).resolve().parent / 'www' / 'js' / 'helper_functions.js'

app_ui = ui.page_fluid(
    # include css and js
    ui.include_css(css_file),
    ui.include_js(js_file),

    # Custom sidebar
    ui.div(id="overlay-bg", onclick="toggleSidebar()"),
    ui.div(
ui.h6(
    ui.img(src='img/settings.svg', height='20px'),
    ' SETTINGS', style="margin-top: 15px;"),
        ui.hr(),
        ui.h6('Statistic settings (entropy, similarity etc)'),
        ui.layout_columns(
            ui.input_numeric('rolling_avg', 'Rolling average', value=1, min=1),
            ui.input_selectize('stat_color', 'Color', list(matplotlib.colors.cnames.keys()),
                               selected='grey'),
        ),
        ui.hr(),
        ui.h6('Alignment settings'),
        ui.row(
            ui.column(
                4,
                ui.input_switch('show_gaps', 'Gaps', value=True),
                ui.input_switch('fancy_gaps', 'Fancy gaps', value=False),
                ui.input_switch('show_legend', 'Legend', value=True),
                ui.input_switch('show_mask', 'Show mask', value=True),
                ui.input_switch('show_ambiguities', 'Show ambiguities', value=True),

            ),
            ui.column(
                4,
                ui.input_selectize('reference', 'Reference', ['first', 'consensus'], selected='first'),
                ui.input_selectize('reference_color', label='Reference color',
                                   choices=list(matplotlib.colors.cnames.keys()), selected='lightsteelblue'),
                ui.input_switch('show_sequence', 'show sequence', value=True),
            ),
            ui.column(
            4,

                ui.input_selectize('matrix', 'Matrix', ['None']),

                ui.input_selectize('matrix_color_mapping', 'Colormap', choices=list(colormaps.keys()),
                                   selected='PuBu_r'),


                ui.input_switch('seq_names', 'show names', value=False),
            )
        ),
        ui.hr(),
        ui.h6('Annotation settings'),
        ui.row(
            ui.column(
            4,
            ui.h6('SNP plot'),
                ui.input_numeric('head_size', 'Head size', value=3, min=1),
                ui.input_numeric('stem_size', 'Stem length', value=1, min=1),
                ui.input_switch('show_legend_third_plot', 'Legend', value=True),
            ),
            ui.column(
                4,
                ui.h6('ORF plot'),
                ui.input_numeric('min_orf_length', 'Length', value=150, min=1),
                ui.input_selectize('color_mapping', 'Colormap', choices=list(colormaps.keys()), selected='jet'),
                ui.input_switch('non_overlapping', 'non-overlapping', value=False),
            ),
            ui.column(
            4,
            ui.h6('Annotation plot'),
                ui.input_selectize('feature_display', 'Feature', ['None']),
                ui.input_selectize('feature_color', 'Color', list(matplotlib.colors.cnames.keys()),
                                   selected='grey'),
                ui.input_numeric('strand_marker_size', 'Strand marker size', value=5, min=1, max=20)
            ),
        ),
        id='overlay-sidebar',
    ),
    # Main side
    ui.navset_bar(
        ui.nav_panel(
            ' UPLOAD/DOWNLOAD',
            ui.div(
                ui.layout_columns(
                    ui.card(
                        ui.card_header(ui.h6('Upload files:')),
                        ui.layout_columns(
                            ui.input_file('alignment_file', 'Multiple sequence alignment:', multiple=False, accept=['.fa', '.fasta', '.aln']),
                            ui.input_file('annotation_file', 'Optional annotation file:', multiple=False, accept=['.gff', '.gff3', '.bed', '.gb']),
                        )
                    ),
                    ui.card(
                        ui.card_header(ui.h6('Download files:'),
                            ui.popover(
                                ui.span(
                                    ui.HTML('<img src="img/gear.svg" alt="settings" style="height:16px; width:16px; position:absolute; top: 10px; right: 7px;">')
                                ),
                                ui.input_selectize('download_type_options', label='Additional options:', choices=['None']),
                                ui.input_selectize('download_format', label='Format:', choices=[]),
                            )
                        ),
                        ui.input_selectize('download_type', label='Choose:', choices=['SNPs']),
                        ui.download_button(
                            'download_stats',
                            'Download',
                            icon=ui.HTML(
                                '<img src="img/download.svg" alt="download icon" style="height:16px; width:16px;">')
                        )
                    )
                ),
                style="display: flex; flex-direction: column; justify-content: center; align-items: center",
            ),
            ui.card(
                ui.h6('About MSAexplorer:'),
                ui.p(
                    "MSAexplorer is an interactive visualization tool designed for exploring multiple sequence alignments (MSAs)."),
                ui.p(ui.a("🔗 Learn more and contribute on GitHub", href="https://github.com/jonas-fuchs/MSAexplorer",
                          target="_blank")),
                class_="about-card"
            ),
            icon=ui.HTML('<img src="img/upload.svg" alt="Upload Icon" style="height: 1em; vertical-align: middle">')
        ),
        ui.nav_panel(
            ' PLOT',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_slider('increase_height', 'Plot height', min=0.5, max=10, step=0.5, value=1),
                    ui.input_selectize('stat_type', ui.h6('First plot'), ['Off'], selected='Off'),
                    ui.input_numeric('plot_1_size', 'Plot fraction',1, min=1, max=200),
                    ui.input_selectize( 'alignment_type', ui.h6('Second plot'), ['Off', 'identity', 'colored identity', 'similarity'], selected='identity'),
                    ui.input_numeric('plot_2_size', 'Plot fraction', 1, min=1, max=200),
                    ui.input_selectize('annotation', ui.h6('Third plot'), ['Off'], selected='Off'),
                    ui.input_numeric('plot_3_size', 'Plot fraction', 1, min=1, max=200),
                    ui.download_button(
                        'download_pdf',
                        'PDF',
                        icon=ui.HTML('<img src="img/download.svg" alt="download icon" style="height:16px; width:16px;">')
                    ),
                    title=ui.h6('Plotting layout'),
                ),
            ui.input_slider('zoom_range', ui.h6('Zoom'), min=0, max=1000, value=(0, 1000), step=1, width='100vw', ticks=True),
                ui.output_plot('msa_plot', height='100vh', width='92vw'),
                fillable=False
            ),
            icon=ui.HTML('<img src="img/chart.svg" alt="Chart Icon" style="height: 1em; vertical-align: middle;">')
        ),
        ui.nav_panel(
        ' ANALYSIS',
            ui.layout_columns(
        ui.value_box(
                    'Alignment type:',
                    ui.output_ui('aln_type'),
                    showcase=ui.HTML('<img src="img/question.svg" style="height:2.5rem; width:2.5rem">'),
                    theme=ui.value_box_theme(bg='#f1f1f1')
                ),
                ui.value_box(
                    'Position range:',
                    ui.output_ui('zoom_range_analysis'),
                    showcase=ui.HTML('<img src="img/arrow_range.svg" style="height:3rem; width:3rem">'),
                ),
                ui.value_box(
                    'Alignment length:',
                    ui.output_ui('aln_len'),
                    showcase=ui.HTML('<img src="img/ruler.svg" style="height:3.5rem; width:3.5rem">'),
                    theme=ui.value_box_theme(bg='#f1f1f1')
                ),
                ui.value_box(
                    'N° of sequences:',
                    ui.output_ui("number_of_seq"),
                    showcase=ui.HTML('<img src="img/number.svg" style="height:2.5rem; width:2.5rem">'),
                ),
                ui.value_box(
                    'Percentage of gaps:',
                    ui.output_ui("per_gaps"),
                    showcase=ui.HTML('<img src="img/percent.svg" style="height:height:2rem; width:2rem">'),
                    theme=ui.value_box_theme(bg='#f1f1f1')
                ),
                ui.value_box(
                    'N° of positions with SNPs:',
                    ui.output_ui("snps"),
                    showcase=ui.HTML('<img src="img/number.svg" style="height:height:2.5rem; width:2.5rem">'),
                ),
            ),
        ui.row(
            ui.column(
                2,
                ui.input_selectize('analysis_plot_type', ui.h6('Left plot'), ['Off', 'Pairwise identity'], selected='Off'),
                ui.input_selectize('additional_analysis_options', ui.h6('Options'), ['None'], selected='None'),
                ui.output_text_verbatim('analysis_info', placeholder=False)
            ),
            ui.column(
                5,
                output_widget('analysis_custom_heatmap'),
            ),
            ui.column(
                5,
                output_widget('analysis_char_freq_heatmap', fillable=True, fill=True),
            ),
        ),
        icon=ui.HTML('<img src="img/analyse.svg" alt="Chart Icon" style="height: 1em; vertical-align: middle;">')
        ),
        ui.nav_spacer(),
        ui.nav_control(
            ui.input_action_button(
            "open_sidebar", "SETTINGS", onclick="toggleSidebar()",
                icon=ui.HTML('<img src="img/settings.svg" alt="Setting Icon" style="height: 1em; vertical-align: middle">'),
            )
        ),
        title=ui.a(
            ui.img(src='img/logo.svg', height='60px'),
            href='https://github.com/jonas-fuchs/MSAexplorer'
        )
    )
)

# set the alignment with minimal user inputs
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


# Define the plotting function
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
    if inputs['stat_type'] != 'Off':
        height_ratios.append(inputs['plot_1_size'])
        plot_functions.append(
            lambda ax: draw.stat_plot(
                aln, ax,
                stat_type=inputs['stat_type'],
                line_width=1,
                rolling_average=inputs['rolling_average'],
                show_x_label=True if inputs['annotation'] == 'Off' and inputs['alignment_type'] == 'Off' else False,
                line_color=inputs['stat_color'])
        )

    # Second plot
    if inputs['alignment_type'] != 'Off':
        height_ratios.append(inputs['plot_2_size'])
        plot_functions.append(
            lambda ax: draw.identity_alignment(
                aln, ax,
                show_sequence=inputs['show_sequence'],
                fancy_gaps=inputs['fancy_gaps'],
                show_gaps=inputs['show_gaps'],
                show_mask=inputs['show_mask'],
                show_mismatches=True,
                show_ambiguities=inputs['show_ambiguities'],
                color_mismatching_chars=True if inputs['alignment_type'] == 'colored identity' else False,
                reference_color=inputs['reference_color'],
                show_seq_names=inputs['seq_names'],
                show_x_label=True if inputs['annotation'] == 'Off' else False,
                show_legend=inputs['show_legend']
        ) if inputs['alignment_type'] == 'identity' or inputs[
            'alignment_type'] == 'colored identity' else draw.similarity_alignment(
                aln, ax,
                show_sequence=inputs['show_sequence'],
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
                aln, ann, ax,
                feature_to_plot=inputs['feature_display'],
                color=inputs['feature_color'],
                direction_marker_size=inputs['strand_marker_size'],
                show_x_label=True
        ) if inputs['annotation'] == 'Annotation' and inputs['annotation_file'] else draw.orf_plot(
                aln, ax,
                cmap=inputs['color_mapping'],
                non_overlapping_orfs=inputs['non_overlapping'],
                direction_marker_size=inputs['strand_marker_size'],
                show_x_label=True,
                show_cbar=inputs['show_legend_third_plot'],
                cbar_fraction=0.2,
                min_length=inputs['min_orf_length']
        ) if inputs['annotation'] == 'Conserved ORFs' else draw.variant_plot(
                aln, ax,
                show_x_label=True,
                lollisize=(inputs['stem_size'], inputs['head_size']),
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


# Define the analysis plotting function (interactive)
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

    if inputs['analysis_plot_type'] == 'Pairwise identity' and inputs['additional_analysis_options'] != 'None':
        matrix = aln.calc_pairwise_identity_matrix(inputs['additional_analysis_options'])
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
                colorbar=dict(thickness=20, ticklen=4, len=0.9, title='%')
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
    Create a frequency heatmap (right heatmap)
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
                colorbar=dict(thickness=20, ticklen=4, title='%')
            )
        )

    fig.update_layout(
        title='Character frequencies (% ungapped)',
        height=figure_height - 50,
        margin=dict(t=50, b=50, l=50, r=50),
    )

    return fig

# Define the server logic
def server(input, output, session):

    reactive.alignment = reactive.Value(None)
    reactive.annotation = reactive.Value(None)

    # create inputs for plotting and pdf
    def prepare_inputs():
        """
        Collect inputs from the UI for the plot tab, only adding the ones needed for enabled features.
        This ensures unnecessary reactivity.
        """
        inputs = {}
        aln = reactive.alignment.get()

        if aln is None:
            return None

        # inhibit the window dimensions and height to trigger a re-rendering
        # as this is only need for the calc whether to show the
        # sequence but itself does not change the plots appearance
        # --> other reactive values will change anyway
        with reactive.isolate():
            dims = input.window_dimensions()
            window_width = dims['width']
            window_height = dims['height']
            increase_height = input.increase_height()

        # Always needed inputs
        inputs['zoom_range'] = input.zoom_range()
        inputs['reference'] = input.reference()
        inputs['alignment_type'] = input.alignment_type()
        inputs['annotation'] = input.annotation()
        inputs['stat_type'] = input.stat_type()

        # STATISTICS (first plot)
        if inputs['stat_type'] != 'Off':
            inputs['plot_1_size'] = input.plot_1_size()
            inputs['rolling_average'] = input.rolling_avg()
            inputs['stat_color'] = input.stat_color()

        # ALIGNMENT (second plot)
        if inputs['alignment_type'] != 'Off':
            inputs['reference_color'] = input.reference_color()
            inputs['plot_2_size'] = input.plot_2_size()
            inputs['fancy_gaps'] = input.fancy_gaps()
            inputs['show_mask'] = input.show_mask()
            inputs['show_gaps'] = input.show_gaps()
            inputs['show_legend'] = input.show_legend()
            inputs['show_ambiguities'] = input.show_ambiguities()
            if inputs['alignment_type'] not in ['identity', 'colored identity']:
                inputs['matrix'] = input.matrix()
                inputs['matrix_color_mapping'] = input.matrix_color_mapping()

            # determine if it makes sense to show the sequence or sequence names
            # therefore figure out if there are enough chars/size that sequence fits in there
            complete_size = input.plot_2_size()
            if inputs['stat_type'] != 'Off':
                complete_size += input.plot_1_size()
            if inputs['annotation'] != 'Off':
                complete_size += input.plot_3_size()
            relative_msa_height = input.plot_2_size() * increase_height / complete_size * window_height / len(aln.alignment)
            relative_msa_width = window_width / (inputs['zoom_range'][1]-inputs['zoom_range'][0])
            # and then decide how to set the show sequence input
            if relative_msa_width >= 11 and relative_msa_height >= 18:
                inputs['show_sequence'] = input.show_sequence()
            else:
                inputs['show_sequence'] = False
            # and the seq_names input - check if there is enough y space for the text
            if relative_msa_height >= 15:
                inputs['seq_names'] = input.seq_names()
            else:
                inputs['seq_names'] = False

        # ANNOTATION (third plot)
        if inputs['annotation'] != 'Off':
            inputs['plot_3_size'] = input.plot_3_size()
            if inputs['annotation'] == 'Annotation':
                inputs['annotation_file'] = input.annotation_file()
                inputs['feature_display'] = input.feature_display()
                inputs['feature_color'] = input.feature_color()
                inputs['strand_marker_size'] = input.strand_marker_size()
            elif inputs['annotation'] == 'Conserved ORFs':
                inputs['color_mapping'] = input.color_mapping()
                inputs['non_overlapping'] = input.non_overlapping()
                inputs['min_orf_length'] = input.min_orf_length()
                inputs['strand_marker_size'] = input.strand_marker_size()
                inputs['show_legend_third_plot'] = input.show_legend_third_plot()
            else:
                inputs['stem_size'] = input.stem_size()
                inputs['head_size'] = input.head_size()
                inputs['show_legend_third_plot'] = input.show_legend_third_plot()

        return inputs

    # separate function if not all inputs are needed
    def prepare_minimal_inputs(zoom: bool = True, window_size: bool = False, ref: bool = False, plot: bool = False):
        """
        minimal inputs with options depending on which are needed
        """
        inputs = {}
        aln = reactive.alignment.get()

        if aln is None:
            return None

        if zoom:
            inputs['zoom_range'] = input.zoom_range()
        if ref:
            inputs['reference'] = input.reference()
        if plot:
            inputs['analysis_plot_type'] = input.analysis_plot_type()
            if inputs['analysis_plot_type'] != 'Off':
                inputs['additional_analysis_options'] = input.additional_analysis_options()
        if window_size:
            inputs['dimensions'] = input.window_dimensions()

        return inputs

    def read_in_annotation(annotation_file):
        """
        Read in an annotation and update the ui accordingly
        """
        ann = explore.Annotation(reactive.alignment.get(), annotation_file[0]['datapath'])
        reactive.annotation.set(ann)

        # update features to display
        ui.update_selectize(
            id='feature_display',
            choices=list(ann.features.keys()),
            selected=list(ann.features.keys())[0]
        )

        # update possible user inputs
        if reactive.alignment.get().aln_type == 'AA':
            ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Annotation'], selected='Annotation')
        else:
            ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Conserved ORFs', 'Annotation'],
                                selected='Annotation')

    # Updates the plot container
    @reactive.Effect
    async def update_height():
        """
        update the plot container height -> Sends a message that is picked up by the js and updates the CSS height
        property for the plot container
        """
        new_plot_height = f'{input.increase_height() * 100}vh'

        await session.send_custom_message("update-plot-container-height", {'height': new_plot_height})

    # Inputs
    @reactive.Effect
    @reactive.event(input.alignment_file)
    def load_alignment():
        """
        Load an alignment and update different now accessible ui elements
        """
        try:
            alignment_file = input.alignment_file()
            annotation_file = input.annotation_file()
            if alignment_file:
                aln = explore.MSA(alignment_file[0]['datapath'], reference_id=None, zoom_range=None)

                # set standard ref
                aln.reference_id = list(aln.alignment.keys())[0]
                reactive.alignment.set(aln)

                # Update zoom slider based on alignment length and user input
                alignment_length = len(next(iter(aln.alignment.values())))-1
                ui.update_slider('zoom_range', max=alignment_length-1, value=(0, alignment_length -1))

                # Update reference
                ui.update_selectize(
                    id='reference', choices=['first', 'consensus'] + list(aln.alignment.keys()), selected='first'
                )

                # Update substitution matrix
                ui.update_selectize(
                    id='matrix',
                    choices=list(config.SUBS_MATRICES[aln.aln_type].keys()),
                    selected='BLOSUM65' if aln.aln_type == 'AA' else 'TRANS',
                )

                # update plot size sliders
                # Adjust the size depending on the number of alignment sequences
                aln_len, seq_threshold = len(aln.alignment.keys()), 5
                for ratio in config.STANDARD_HEIGHT_RATIOS.keys():
                    if aln_len >= ratio:
                        seq_threshold = ratio

                # update standard settings
                ui.update_numeric('plot_1_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][0])
                ui.update_numeric('plot_2_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][1])
                ui.update_numeric('plot_3_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][2])
                # Some of the function highly depend on the alignment type
                if aln.aln_type == 'AA':
                    ui.update_selectize('stat_type', choices=['Off', 'entropy', 'coverage', 'identity', 'similarity'], selected='Off')
                    ui.update_selectize('annotation', choices=['Off', 'SNPs'])
                else:
                    # needed because if an as aln and then a nt aln are loaded it will not change
                    ui.update_selectize('stat_type', choices=['Off', 'gc', 'entropy', 'coverage', 'identity', 'similarity', 'ts tv score'], selected='Off')
                    ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Conserved ORFs'])
                # case if annotation file is uploaded prior to the alignment file
                if annotation_file:
                    read_in_annotation(annotation_file)
        # show the user if something with parsing went wrong
        except Exception as e:
            print(f"Error: {e}")  # print to console
            ui.notification_show(ui.tags.div(
                    f'Error: {e}',
                    style="color: red; font-weight: bold;"
                ), duration=10)  # print to user

    @reactive.Effect
    @reactive.event(input.annotation_file)
    def load_annotation():
        """
        load the annotation - catches errors if its the wrong format and displays it (otherwise app would crash)
        """
        # read annotation file
        try:
            annotation_file, alignment_file = input.annotation_file(), input.alignment_file()
            if annotation_file and alignment_file:
                read_in_annotation(annotation_file)
        # show the user if something with parsing went wrong
        except Exception as e:
            print(f"Error: {e}")  # print to console
            ui.notification_show(ui.tags.div(
                f'Error: {e}',
                style="color: red; font-weight: bold;"
            ), duration=10)  # print to user


    @reactive.Effect
    @reactive.event(input.download_type)
    def update_download_options():
        """
        Update the UI - which download options are displayed
        """
        if input.download_type() == 'SNPs':
            ui.update_selectize('download_format', choices=['vcf', 'tabular'])

    @render.download()
    def download_stats():
        """
        Download various files in standard format
        """
        try:
            # Initialize
            download_format = input.download_format()
            aln = reactive.alignment.get()
            if aln is None:
                raise FileNotFoundError("No alignment data available. Please upload an alignment.")

            if input.download_type() == 'SNPs':
                download_data = export.snps(aln.get_snps(), format_type=download_format)
                prefix = 'SNPs_'

            # Create a temporary file for the download
            with tempfile.NamedTemporaryFile(prefix=prefix, suffix=f'.{download_format}', delete=False) as tmpfile:
                tmpfile.write(download_data.encode('utf-8'))
                tmpfile.flush()  # Ensure data is written to disk

                return tmpfile.name
        # also send a notification to the user
        except FileNotFoundError:
            ui.notification_show(ui.tags.div(
                'No alignment was uploaded.',
                style="color: red; font-weight: bold;"
            ), duration=10)


    # Outputs
    @output
    @render.plot
    def msa_plot():
        """
        plot the alignment
        """
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        return create_msa_plot(aln, ann, prepare_inputs())

    @output
    @render.download
    def download_pdf():
        """
        allows the export of the msa plot as pdf - importantly it will be same size as the
        window ensuring that the plot is exactly the same as the rendered displayed plot
        """
        # get annotation and alignment
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        # access the window dimensions
        dimensions = input.window_dimensions_plot()
        figure_width_inches = dimensions['width'] / 96
        figure_height_inches = dimensions['height'] / 96 * input.increase_height()

        # plot with a temp name
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            fig = create_msa_plot(aln, ann, prepare_inputs(), fig_size=(figure_width_inches, figure_height_inches))
            # tight layout needed here to plot everything correctly
            fig.tight_layout()
            fig.savefig(tmpfile.name, format="pdf")
            plt.close(fig)
            return tmpfile.name

    # showcases:
    @render.ui
    def aln_type():
        aln = reactive.alignment.get()
        if aln is None:
            return None
        return aln.aln_type

    @render.ui
    def zoom_range_analysis():
        aln = reactive.alignment.get()
        if aln is None:
            return None

        aln = set_aln(aln, prepare_minimal_inputs())

        return f'{aln.zoom[0]} - {aln.zoom[1]}'

    @render.ui
    def number_of_seq():
        aln = reactive.alignment.get()
        if aln is None:
            return None

        return len(aln.alignment)

    @render.ui
    def aln_len():
        aln = reactive.alignment.get()
        if aln is None:
            return None

        aln = set_aln(aln, prepare_minimal_inputs())

        return aln.length

    @render.ui
    def per_gaps():
        aln = reactive.alignment.get()
        if aln is None:
            return None

        aln = set_aln(aln, prepare_minimal_inputs())

        return round(aln.calc_character_frequencies()['total']['-']['% of alignment'], 2)

    @render.ui
    def snps():
        aln = reactive.alignment.get()
        if aln is None:
            return None

        aln = set_aln(aln, prepare_minimal_inputs(ref=True))

        return len(aln.get_snps()['POS'])

    @reactive.Effect
    @reactive.event(input.analysis_plot_type)
    def update_additional_options():
        """
        Update UI for the analysis tab
        """
        # ensure that it is switched back
        if input.analysis_plot_type() == 'Off':
            ui.update_selectize('additional_analysis_options', choices=['None'], selected='None')
        if input.analysis_plot_type() == 'Pairwise identity':
            ui.update_selectize('additional_analysis_options', choices=['ghd', 'lhd', 'ged', 'gcd'], selected='ghd')

    @render.text
    def analysis_info():
        """
        show custom text for additional options
        """
        selected_option = input.additional_analysis_options()
        if selected_option == "ghd":
            return 'INFO ghd (global hamming distance):\n\nAt each alignment position, check if\ncharacters match:\n\ndistance = matches / alignment_length * 100'
        elif selected_option == "lhd":
            return 'INFO lhd (local hamming distance):\n\nRestrict the alignment to the region\nin both sequences that do not start\nand end with gaps:\n\ndistance = matches / min(end-ungapped seq1, end-ungapped seq2) * 100'
        elif selected_option == 'ged':
            return 'INFO ged (gap excluded distance):\n\nAll gaps are excluded from the \nalignment\n\ndistance = matches / (matches + mismatches) * 100'
        elif selected_option == 'gcd':
            return 'INFO gcd (gap compressed distance):\n\nAll consecutive gaps arecompressed to\none mismatch.\n\ndistance = matches / gap_compressed_alignment_length * 100'
        else:
            return None


    @render_widget
    def analysis_custom_heatmap():
        """
        Create the heatmap
        """
        aln = reactive.alignment.get()

        return create_analysis_custom_heatmap(aln, prepare_minimal_inputs( plot=True, window_size=True))


    @render_widget
    def analysis_char_freq_heatmap():
        aln = reactive.alignment.get()

        return create_freq_heatmap(aln,  prepare_minimal_inputs(window_size=True))


# run the app
app = App(app_ui, server, static_assets={'/img': Path(__file__).parent/'img'})
