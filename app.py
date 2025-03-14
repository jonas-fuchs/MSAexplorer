"""
This contains the code to create the MSAexplorer shiny application
"""

# build-in
from pathlib import Path
import tempfile
from functools import lru_cache

# libs
from shiny import App, render, ui, reactive
import matplotlib
from matplotlib import colormaps
import matplotlib.pyplot as plt

# msaexplorer
from msaexplorer import explore, draw, config

# file paths for css and js
css_file = Path(__file__).resolve().parent / 'www' /'css' / 'styles.css'
js_file = Path(__file__).resolve().parent / 'www' / 'js' / 'window_dimensions.js'

# define the UI
app_ui = ui.page_fluid(
    ui.img(
            src='assets/logo.svg', height='100px'
    ),
    # include css
    ui.include_css(css_file),
    # get the input dimensions (separate js)
    ui.include_js(js_file),
    ui.div(
        ui.a(
            ui.img(src='https://github.githubassets.com/assets/GitHub-Logo-ee398b662d42.png', height='15px'),
            href='https://github.com/jonas-fuchs/MSAexplorer',
            target='_blank',
            title='Give it a star!'
        ),
        style='position: fixed; top: 10px; right: 10px; z-index: 1000;'
    ),
    ui.navset_tab(
        ui.nav_panel(
            'Upload Files',
            ui.tooltip(
                ui.input_file('alignment_file', ui.h6('MSA', class_='section-title'), multiple=False),
                'Multiple sequence alignment file to display.'
            ),
            ui.tooltip(
                ui.input_file('annotation_file', ui.h6('Optional .gff3, .bed or .gb', class_='section-title'), multiple=False),
                'Optional annotation file to display. Sequence id must be present in the alignment for correct mapping.'
            )
        ),
        ui.nav_panel(
        'Advanced Settings',
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
                        'Whether to show gaps as an actual gap.'
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
                    )
                ),
            )
        ),
        ui.nav_panel(
            'Visualization',
            ui.layout_sidebar(
                ui.sidebar(
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
                    ),
                    ui.tooltip(
                        ui.input_switch('show_sequence', 'show sequence', value=False),
                        'Whether to show the sequence if zoomed.'
                    ),
                    ui.tooltip(
                        ui.input_switch('seq_names', 'show names', value=False),
                        'Whether to show sequence names at the left side of the alignment'
                    ),
                    ui.tooltip(
                        ui.download_button('download_pdf', 'PDF'),
                        'Get the plot as a pdf.'
                    )
                ),
                ui.output_plot('msa_plot', height='100vh', width='92vw'),
                ui.input_slider('zoom_range', ui.h6('Zoom'), min=0, max=1000, value=(0, 1000), step=1, width='100vw', ticks=True),
            ),
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
                show_sequence=inputs['show_sequence'] if aln.length/inputs['window_width'] <= 0.085 else False,
                fancy_gaps=inputs['show_gaps'],
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
                show_sequence=inputs['show_sequence'] if aln.length/inputs['window_width'] <= 0.085 else False,
                fancy_gaps=inputs['show_gaps'],
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

    # create inputs for plotting and pdf
    def prepare_inputs():
        # Collect inputs from the UI
        return {
            'reference': input.reference(),
            'reference_color': input.reference_color(),
            'show_mask': input.show_mask(),
            'show_ambiguities': input.show_ambiguities(),
            'window_width': input.window_dimensions()['width'],
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
            'seq_names': input.seq_names(),
            'show_sequence': input.show_sequence(),
            'show_legend': input.show_legend(),
            'annotation': input.annotation(),
            'annotation_file': input.annotation_file(),
            'feature_display': input.feature_display(),
            'color_mapping': input.color_mapping(),
            'non_overlapping': input.non_overlapping(),
            'min_orf_length': input.min_orf_length(),
            'stem_size': input.stem_size(),
            'head_size': input.head_size(),
            'show_legend_third_plot': input.show_legend_third_plot(),
        }


    @reactive.Effect
    @reactive.event(input.alignment_file)
    def load_alignment():
        alignment_file = input.alignment_file()
        if alignment_file:
            aln = explore.MSA(alignment_file[0]['datapath'], reference_id=None, zoom_range=None)

            # set standard ref
            aln.reference_id = list(aln.alignment.keys())[0]
            reactive.alignment.set(aln)

            # Update zoom slider based on alignment length and user input
            alignment_length = len(next(iter(aln.alignment.values())))-1
            ui.update_slider('zoom_range', max=alignment_length-1, value=(0, alignment_length-1))

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


    @reactive.Effect
    @reactive.event(input.annotation_file)
    def load_annotation():
        # read annotation file
        annotation_file = input.annotation_file()
        if annotation_file:
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
                ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Conserved ORFs', 'Annotation'], selected='Annotation')

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
        figure_height_inches = dimensions['height'] / 96

        # plot with a temp name
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            fig = create_msa_plot(aln, ann, prepare_inputs(), fig_size=(figure_width_inches, figure_height_inches))
            # tight layout needed here to plot everything correctly
            fig.tight_layout()
            fig.savefig(tmpfile.name, format="pdf")
            plt.close(fig)
            return tmpfile.name

# run the app
app = App(app_ui, server, static_assets={'/assets': Path(__file__).parent/"assets"})


#TODO: Analysis tab