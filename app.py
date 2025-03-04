"""
This contains the code to create the MSAexplorer shiny application
"""

# build-in
from pathlib import Path
import tempfile

# libs
from shiny import App, render, ui, reactive
import matplotlib
from matplotlib import colormaps
import matplotlib.pyplot as plt

# msaexplorer
from msaexplorer import explore, draw, config

# file paths for css and js
css_file = Path(__file__).parent / 'www' /'css' / 'styles.css'
js_file = Path(__file__).parent / 'www' / 'js' / 'window_dimensions.js'

# define the UI
app_ui = ui.page_fluid(
    ui.h2('MSAexplorer'),
    # include css
    ui.include_css(css_file),
    # get the input dimensions (separate js)
    ui.include_js(js_file),
    ui.div(
        ui.a(
            ui.img(src='https://github.githubassets.com/assets/GitHub-Logo-ee398b662d42.png', height='15px'),
            href='https://github.com/jonas-fuchs/MSAexplorer',
            target='_blank',
            title='View on GitHub'
        ),
        style='position: fixed; top: 10px; right: 10px; z-index: 1000;'
    ),
    ui.navset_tab(
        ui.nav_panel(
            'Upload Files',
            ui.input_file('alignment_file', ui.h6('Multiple sequence alignment', class_='section-title'), multiple=False),
            ui.input_file('annotation_file', ui.h6('Optional .gff3, .bed or .gb', class_='section-title'), multiple=False),
        ),
        ui.nav_panel(
        'Advanced Settings',
            ui.row(
                ui.h6('First plot', class_='section-title'),
                ui.column(
                    4,
                    ui.input_numeric('rolling_avg', 'Rolling average', value=20, min=1)
                ),
                ui.column(
                    4,
                    ui.input_selectize('stat_color', 'Line color', list(matplotlib.colors.cnames.keys()), selected='grey')
                )
            ),
            ui.row(
                ui.h6('Second plot', class_='section-title'),
                ui.column(
                    4,
                    ui.input_switch('show_gaps', 'Show gaps', value=True),
                    ui.input_switch('show_legend', 'Show legend', value=True)
                ),
                ui.column(
                    4,
                    ui.input_selectize('reference', 'Reference sequence', ['first' ,'consensus'], selected='first'),
                ),
                ui.column(
                    4,
                    ui.input_selectize('matrix', 'Substitution matrix (similarity)', ['None'])
                )
            ),
            ui.row(
                ui.h6('Third plot',  class_='section-title'),
                ui.column(
                    4,
                ui.h6('SNP plot'),
                    ui.input_numeric('head_size', 'Variant size (Head)', value=3, min=1),
                    ui.input_numeric('stem_size', 'Variant size (Stem)', value=1, min=1),
                    ui.input_switch('show_legend_variants', 'Show legend', value=True)
                ),
                ui.column(
                    4,
                ui.h6('ORF plot'),
                    ui.input_numeric('min_orf_length', 'Minimum ORF length', value=150, min=1),
                    ui.input_switch('non_overlapping', 'Non-Overlapping ORFs', value=False),
                    ui.input_selectize('color_mapping', 'Colormap for conservation', choices=list(colormaps.keys()), selected='jet'
                                       ),
                ),
                ui.column(
                    4,
                ui.h6('Annotation plot'),
                    ui.input_selectize('feature_display', 'Feature to display', ['None']),
                    ui.input_selectize('feature_color', 'Feature color', list(matplotlib.colors.cnames.keys()), selected='grey')
                ),
            )
        ),
        ui.nav_panel(
            'Visualization',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize('stat_type', ui.h6('First plot'), ['Off', 'gc', 'entropy', 'coverage', 'identity'], selected='gc'),
                    ui.input_numeric('plot_1_size', 'Plot fraction',1, min=1, max=200),
                    ui.input_selectize( 'alignment_type', ui.h6('Second plot'), ['Off', 'identity', 'similarity'], selected='identity'),
                    ui.input_numeric('plot_2_size', 'Plot fraction', 1, min=1,step=2, max=200),
                    ui.input_selectize('annotation', ui.h6('Third plot'), ['Off', 'SNPs','Conserved ORFs', 'Annotation'], selected='Annotation'),
                    ui.input_numeric('plot_3_size', 'Plot fraction', 1, min=1, max=200),
                    ui.input_slider('zoom_range', ui.h6('Zoom'), min=0, max=1000, value=(0, 1000), step=1),
                    ui.input_switch('seq_names', 'show names', value=False),
                    ui.download_button('download_pdf', 'PDF')
                ),
                ui.output_plot('msa_plot', height='90vh', width='90vw'),
            ),
        )
    )
)

# Define the server logic
def server(input, output, session):

    reactive.alignment = reactive.Value(None)
    reactive.annotation = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.alignment_file)
    def load_alignment():
        alignment_file = input.alignment_file()
        if alignment_file:
            aln = explore.MSA(alignment_file[0]['datapath'], reference_id=None, zoom_range=None)
            # set standard ref
            aln.reference_id = list(aln.alignment.keys())[0]
            reactive.alignment.set(aln)
            # Update zoom slider based on alignment length
            alignment_length = len(next(iter(aln.alignment.values())))-1
            # Update zoom slider
            ui.update_slider('zoom_range', max=alignment_length, value=(0, alignment_length))
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

            ui.update_numeric('plot_1_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][0])
            ui.update_numeric('plot_2_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][1])
            ui.update_numeric('plot_3_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][2])


    @reactive.Effect
    @reactive.event(input.annotation_file)
    def load_annotation():
        annotation_file = input.annotation_file()
        if annotation_file:
            ann = explore.Annotation(reactive.alignment.get(), annotation_file[0]['datapath'])
            reactive.annotation.set(ann)  # load alignment
            # update features to display
            ui.update_selectize(
                id='feature_display',
                choices=list(ann.features.keys()),
                selected=list(ann.features.keys())[0]
            )

    def create_msa_plot(aln, ann, inputs, fig_size=None):

        if not aln:
            return None

        # set the reference sequence
        if input.reference() == 'first':
            aln.reference_id = list(aln.alignment.keys())[0]
        elif input.reference() == 'consensus':
            aln.reference_id = None
        else:
            aln.reference_id = inputs['reference']

        # Update zoom level from slider
        aln.zoom = tuple(inputs['zoom_range'])

        # Prepare the plot with 3 subplots
        # Fig_size is needed for pdf download
        if fig_size is not None:
            fig, axes = plt.subplots(nrows=3, height_ratios=[inputs['plot_1_size'], inputs['plot_2_size'], inputs['plot_3_size']], figsize=fig_size)
        else:
            fig, axes = plt.subplots(nrows=3, height_ratios=[inputs['plot_1_size'], inputs['plot_2_size'], inputs['plot_3_size']])
        # Subplot 1: Stats Plot
        draw.stat_plot(
            aln,
            axes[0],
            stat_type=inputs['stat_type'],
            line_width=1,
            rolling_average=inputs['rolling_average'],
            line_color=inputs['stat_color']
        )

        # Subplot 2: Alignment Plot (Identity or Similarity)
        if inputs['alignment_type'] == 'identity':
            draw.identity_alignment(
                aln, axes[1],
                fancy_gaps=inputs['show_gaps'],
                show_gaps=inputs['show_gaps'],
                show_mask=True,
                show_mismatches=True,
                show_ambiguities=True,
                reference_color='lightsteelblue',
                show_seq_names=inputs['seq_names'],
                show_x_label=True if inputs['annotation'] == 'Annotation' and not inputs['annotation_file'] else False,
                show_legend=inputs['show_legend']
            )
        else:
            draw.similarity_alignment(
                aln, axes[1],
                fancy_gaps=inputs['show_gaps'],
                show_gaps=inputs['show_gaps'],
                matrix_type=inputs['matrix'],
                show_seq_names=inputs['seq_names'],
                show_cbar=inputs['show_legend'],
                cbar_fraction=0.02,
                show_x_label=True if inputs['annotation'] == 'Annotation' and not inputs['annotation_file'] else False
            )

        # Subplot 3: Annotation or ORF Plot
        if inputs['annotation'] == 'Annotation' and inputs['annotation_file']:
            draw.annotation_plot(
                aln, ann,
                axes[2],
                feature_to_plot=inputs['feature_display'],
                show_x_label=True
            )
        elif inputs['annotation'] == 'Conserved ORFs':
            draw.orf_plot(
                aln, axes[2],
                cmap=inputs['color_mapping'],
                non_overlapping_orfs=inputs['non_overlapping'],
                show_x_label=True,
                show_cbar=True,
                cbar_fraction=0.2,
                min_length=inputs['min_orf_length']
            )
        elif inputs['annotation'] == 'SNPs':
            draw.variant_plot(
                aln, axes[2],
                show_x_label=True,
                lollisize=(inputs['stem_size'], inputs['head_size']),
                show_legend=inputs['show_legend_variants']
            )
        # turn everything off
        else:
            axes[2].axis('off')

        return fig

    @output
    @render.plot
    def msa_plot():
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        # Collect inputs from the UI
        inputs = {
            'reference': input.reference(),
            'zoom_range': input.zoom_range(),
            'plot_1_size': input.plot_1_size(),
            'plot_2_size': input.plot_2_size(),
            'plot_3_size': input.plot_3_size(),
            'stat_type': input.stat_type(),
            'rolling_average': input.rolling_avg(),
            'stat_color': input.stat_color(),
            'alignment_type': input.alignment_type(),
            'matrix': input.matrix(),
            'show_gaps': input.show_gaps(),
            'seq_names': input.seq_names(),
            'show_legend': input.show_legend(),
            'annotation': input.annotation(),
            'annotation_file': input.annotation_file(),
            'feature_display': input.feature_display(),
            'color_mapping': input.color_mapping(),
            'non_overlapping': input.non_overlapping(),
            'min_orf_length': input.min_orf_length(),
            'stem_size': input.stem_size(),
            'head_size': input.head_size(),
            'show_legend_variants': input.show_legend_variants(),
        }

        fig = create_msa_plot(aln, ann, inputs)
        return fig

    @output
    @render.download
    def download_pdf():
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        # Access the window dimensions
        dimensions = input.window_dimensions()
        if not dimensions:
            raise ValueError("Window dimensions not available.")

        # Convert window dimensions (pixels) to inches
        screen_dpi = 96  # Typical screen DPI
        figure_width_inches = dimensions['width'] / screen_dpi
        figure_height_inches = dimensions['height'] / screen_dpi


        # Collect inputs from the UI
        inputs = {
            'reference': input.reference(),
            'zoom_range': input.zoom_range(),
            'plot_1_size': input.plot_1_size(),
            'plot_2_size': input.plot_2_size(),
            'plot_3_size': input.plot_3_size(),
            'stat_type': input.stat_type(),
            'rolling_average': input.rolling_avg(),
            'stat_color': input.stat_color(),
            'alignment_type': input.alignment_type(),
            'matrix': input.matrix(),
            'show_gaps': input.show_gaps(),
            'seq_names': input.seq_names(),
            'show_legend': input.show_legend(),
            'annotation': input.annotation(),
            'annotation_file': input.annotation_file(),
            'feature_display': input.feature_display(),
            'color_mapping': input.color_mapping(),
            'non_overlapping': input.non_overlapping(),
            'min_orf_length': input.min_orf_length(),
            'stem_size': input.stem_size(),
            'head_size': input.head_size(),
            'show_legend_variants': input.show_legend_variants(),
        }
        # plot with a temp name
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            fig = create_msa_plot(aln, ann, inputs, fig_size=(figure_width_inches, figure_height_inches))
            # tight layout needed here to plot everything correctly
            fig.tight_layout()
            fig.savefig(tmpfile.name, format="pdf")
            plt.close(fig)
            return tmpfile.name

app = App(app_ui, server)

# TODO: on/off for plots
# TODO: on/off for plots for as alignments

