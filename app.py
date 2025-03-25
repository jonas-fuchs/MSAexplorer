"""
This contains the code to create the MSAexplorer shiny application
"""

# build-in
from pathlib import Path
import tempfile
import sys, os

# libs
from shiny import App, render, ui, reactive
import matplotlib
from matplotlib import colormaps
import matplotlib.pyplot as plt

# msaexplorer
from msaexplorer import explore, draw, config

# file paths for css and js
css_file = Path(__file__).resolve().parent / 'www' /'css' / 'styles.css'
js_file = Path(__file__).resolve().parent / 'www' / 'js' / 'helper_functions.js'

# define the UI
app_ui = ui.page_fluid(
    ui.img(
            src='img/logo.svg', height='100px'
    ),
    # include css
    ui.include_css(css_file),
    # get the input dimensions (separate js)
    ui.include_js(js_file),
    ui.a(
        ui.img(
            src='img/github.svg',
            height='40px',
            class_='github-logo'
        ),
        href='https://github.com/jonas-fuchs/MSAexplorer',
        target='_blank',
        title='Give it a star!',
        style='position: fixed; top: 10px; right: 10px; z-index: 1000;'
    ),
    ui.navset_tab(
        ui.nav_panel(
            ' Upload',
            ui.tooltip(
                ui.input_file('alignment_file', ui.h6('Multiple sequence alignment', class_='section-title'), multiple=False, accept=['.fa', '.fasta', '.aln']),
                'Multiple sequence alignment file to display.'
            ),
            ui.tooltip(
                ui.input_file('annotation_file', ui.h6('Optional annotation file', class_='section-title'), multiple=False, accept=['.gff', '.gff3', '.bed', '.gb']),
                'Optional annotation file to display. Sequence id must be present in the alignment for correct mapping.'
            ),
            icon=ui.HTML('<img src="img/upload.svg" alt="Upload Icon" style="height: 1em; vertical-align: middle">')
        ),
        ui.nav_panel(
        ' Advanced Settings',
            ui.row(
                ui.h6('First plot', class_='section-title'),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.input_numeric('rolling_avg', 'Average', value=1, min=1),
                        'Rolling average over character intervals.'
                    )
                ),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.input_selectize('stat_color', 'Color', list(matplotlib.colors.cnames.keys()), selected='grey'),
                        'Any named matplotlib color for the line.'
                    )
                )
            ),
            ui.row(
                ui.h6('Second plot', class_='section-title'),
                ui.column(
                    4,ui.tooltip(
                        ui.input_switch('show_gaps', 'Gaps', value=True),
                        'Whether to show gaps.'
                    ),
                    ui.tooltip(
                        ui.input_switch('fancy_gaps', 'Fancy gaps', value=False),
                        'Whether to show gaps with a black bar.'
                    ),
                    ui.tooltip(
                        ui.input_switch('show_legend', 'Legend', value=True),
                        'Whether to show legend'
                    ),
                    ui.tooltip(
                        ui.input_switch('show_mask', 'Show mask', value=True),
                        'Show masked charaters "X" or "N" - only relevant for identity alignments.'
                    ),
                    ui.tooltip(
                        ui.input_switch('show_ambiguities', 'Show ambiguities', value=True),
                        'Show masked ambiguities - only relevant for nt identity alignments.'
                    ),
                ),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.input_selectize('reference', 'Reference', ['first', 'consensus'], selected='first'),
                        'Which reference sequence to calculate identity/similarity.'
                    ),
                    ui.tooltip(
                        ui.input_selectize('reference_color', label='Reference color', choices=list(matplotlib.colors.cnames.keys()), selected='lightsteelblue'),
                        'Color for the reference sequence.'
                    )
                ),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.input_selectize('matrix', 'Matrix', ['None']),
                        'Substitution matrix for similarity mapping.'
                    ),
                        ui.tooltip(
                            ui.input_selectize('matrix_color_mapping', 'Colormap', choices=list(colormaps.keys()),
                                               selected='PuBu_r'),
                            'colormap for similarity plots - any matplotlib colormap.'
                    )
                )
            ),
            ui.row(
                ui.h6('Third plot',  class_='section-title'),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.input_switch('show_legend_third_plot', 'Legend', value=True),
                        'Whether to show legend for the third plot.'
                    ),
                )
            ),
            ui.row(
                ui.column(
                    4,
                    ui.h6('SNP plot'),
                    ui.tooltip(
                        ui.input_numeric('head_size', 'Head size', value=3, min=1),
                        'Size of the head dot.'
                    ),
                    ui.tooltip(
                        ui.input_numeric('stem_size', 'Stem length', value=1, min=1),
                        'Length of the stem.'
                    ),
                ),
                ui.column(
                    4,
                    ui.tooltip(
                        ui.h6('ORF plot'),
                        'Only relevant for nt alignments.',
                        placement='left'
                     ),
                    ui.tooltip(
                        ui.input_numeric('min_orf_length', 'Length', value=150, min=1),
                        'Minimum ORF length to calculate.'
                    ),
                    ui.tooltip(
                        ui.input_selectize('color_mapping', 'Colormap', choices=list(colormaps.keys()), selected='jet'),
                        'Colormap for conservation - any matplotlib colormap.'
                    ),
                    ui.tooltip(
                        ui.input_switch('non_overlapping', 'non-overlapping', value=False),
                        'Whether to show non-overlapping ORFs - greedy: works from 5 to 3 prime.'
                    ),
                ),
                ui.column(
                    4,
                ui.h6('Annotation plot'),
                    ui.tooltip(
                        ui.input_selectize('feature_display', 'Feature', ['None']),
                        'Which feature to display.'
                    ),
                    ui.tooltip(
                        ui.input_selectize('feature_color', 'Color', list(matplotlib.colors.cnames.keys()), selected='grey'),
                        'Color of the feature.'
                    ),
                    ui.tooltip(
                        ui.input_numeric('strand_marker_size', 'Strand marker size', value=5, min=1, max=20),
                        'Marker size (triangle) of the strand.'
                    ),
                ),
            ),
            icon=ui.HTML('<img src="img/settings.svg" alt="Setting Icon" style="height: 1em; vertical-align: middle">')
        ),
        ui.nav_panel(
            ' Visualization',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h6('General settings'),
                    ui.tooltip(
                        ui.input_slider(
                            'increase_height', 'Plot height', min=0.5, max=10, step=0.5, value=1,
                        ),
                        'Height relative to your window height.'
                    ),
                    ui.tooltip(
                        ui.input_switch('show_sequence', 'show sequence', value=True),
                        'Whether to show the sequence if zoomed.'
                    ),
                    ui.tooltip(
                        ui.input_switch('seq_names', 'show names', value=False),
                        'Whether to show sequence names at the left side of the alignment'
                    ),
                    ui.tooltip(
                        ui.download_button('download_pdf', 'PDF'),
                        'Get the plot as a pdf.'
                    ),
                    ui.input_selectize('stat_type', ui.h6('First plot'), ['Off'], selected='Off'),
                    ui.tooltip(
                        ui.input_numeric('plot_1_size', 'Plot fraction',1, min=1, max=200),
                        'Fraction of the total plot size'
                    ),
                    ui.input_selectize( 'alignment_type', ui.h6('Second plot'), ['Off', 'identity', 'colored identity', 'similarity'], selected='identity'),
                    ui.tooltip(
                        ui.input_numeric('plot_2_size', 'Plot fraction', 1, min=1, max=200),
                        'Fraction of the total plot size'
                    ),
                    ui.input_selectize('annotation', ui.h6('Third plot'), ['Off'], selected='Off'),
                    ui.tooltip(
                        ui.input_numeric('plot_3_size', 'Plot fraction', 1, min=1, max=200),
                        'Fraction of the total plot size'
                    )
                ),
                ui.output_plot('msa_plot', height='100vh', width='92vw'),
                ui.input_slider('zoom_range', ui.h6('Zoom'), min=0, max=1000, value=(0, 1000), step=1, width='100vw', ticks=True),
            ),
        icon=ui.HTML('<img src="img/chart.svg" alt="Chart Icon" style="height: 1em; vertical-align: middle;">')
        )
    )
)


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

    # set the reference sequence
    if inputs['reference'] == 'first':
        aln.reference_id = list(aln.alignment.keys())[0]
    elif inputs['reference'] == 'consensus':
        aln.reference_id = None
    else:
        aln.reference_id = inputs['reference']

    # Update zoom level from slider -> +1 needed as msaexplorer uses range for zoom
    aln.zoom = (inputs['zoom_range'][0], inputs['zoom_range'][1] + 1)

    # determine height and width where it makes sense to plot the text
    if inputs['alignment_type'] != 'Off':
        complete_size = inputs['plot_2_size']
        if inputs['stat_type'] != 'Off':
            complete_size += inputs['plot_1_size']
        if inputs['annotation'] != 'Off':
            complete_size += inputs['plot_3_size']
        plot_2_height_ratio = inputs['plot_2_size'] * inputs['input_increase_height']/complete_size
        relative_msa_height = plot_2_height_ratio * inputs['window_height'] / len(aln.alignment)
        relative_msa_width = inputs['window_width']/aln.length

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
                show_sequence=inputs['show_sequence'] if relative_msa_width >= 11 and relative_msa_height >= 20 else False,
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
                show_sequence=inputs['show_sequence'] if relative_msa_width >= 11 and relative_msa_height >= 20 else False,
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


# Define the server logic
def server(input, output, session):

    reactive.alignment = reactive.Value(None)
    reactive.annotation = reactive.Value(None)
    height_value = reactive.Value('100vh')
    # create inputs for plotting and pdf
    def prepare_inputs():
        """Collect inputs from the UI"""
        return {
            'reference': input.reference(),
            'reference_color': input.reference_color(),
            'show_mask': input.show_mask(),
            'show_ambiguities': input.show_ambiguities(),
            'window_width': input.window_dimensions()['width'],
            'window_height': input.window_dimensions()['height'],
            'zoom_range': input.zoom_range(),
            'plot_1_size': input.plot_1_size(),
            'plot_2_size': input.plot_2_size(),
            'plot_3_size': input.plot_3_size(),
            'stat_type': input.stat_type(),
            'rolling_average': input.rolling_avg(),
            'stat_color': input.stat_color(),
            'alignment_type': input.alignment_type(),
            'matrix': input.matrix(),
            'matrix_color_mapping': input.matrix_color_mapping(),
            'show_gaps': input.show_gaps(),
            'fancy_gaps': input.fancy_gaps(),
            'seq_names': input.seq_names(),
            'show_sequence': input.show_sequence(),
            'show_legend': input.show_legend(),
            'annotation': input.annotation(),
            'annotation_file': input.annotation_file(),
            'feature_display': input.feature_display(),
            'feature_color': input.feature_color(),
            'strand_marker_size': input.strand_marker_size(),
            'color_mapping': input.color_mapping(),
            'non_overlapping': input.non_overlapping(),
            'min_orf_length': input.min_orf_length(),
            'stem_size': input.stem_size(),
            'head_size': input.head_size(),
            'show_legend_third_plot': input.show_legend_third_plot(),
            'input_increase_height': input.increase_height()
        }

    def read_in_annotation(annotation_file):
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


    @reactive.Effect
    async def update_width():
        new_plot_height = f'{input.increase_height()*100}vh'
        new_sidebar_height = f'{input.increase_height()*(100+40)}vh'
        height_value.set(new_plot_height)

        await session.send_custom_message('update-plot-height', {'height': new_plot_height}),
        await session.send_custom_message('update-sidebar-height', {'height': new_sidebar_height})

    @reactive.Effect
    @reactive.event(input.alignment_file)
    def load_alignment():
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
                ui.update_slider('zoom_range', max=alignment_length-1, value=(0, int(alignment_length/10)))

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

    @output
    @render.plot
    def msa_plot():
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        return create_msa_plot(aln, ann, prepare_inputs())

    @output
    @render.download
    def download_pdf():
        # get annotation and alignment
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        # access the window dimensions
        dimensions = input.window_dimensions()
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

# run the app
app = App(app_ui, server, static_assets={'/img': Path(__file__).parent/'img'})


#TODO: Analysis tab