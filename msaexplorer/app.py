# build in
from pathlib import Path

# libs
from shiny import App, render, ui, reactive
import matplotlib.pyplot as plt

# msaexplorer
from msaexplorer import explore, draw, config

# define the UI
app_ui = ui.page_fluid(
    ui.h2('MSAexplorer'),
    ui.include_css(
        Path(__file__).parent / 'www/styles.css'
    ),
    ui.navset_tab(
        ui.nav_panel(
            'Upload Files',
            ui.input_file('alignment_file', ui.h6('Upload Alignment File (.aln)', class_='section-title'), multiple=False),
            ui.input_file('annotation_file', ui.h6('Annotation File (.gff3, .bed, .gb)', class_='section-title'), multiple=False),
        ),
        ui.nav_panel(
        'Advanced Settings',
            ui.row(
                ui.h6('First plot:', class_='section-title'),
                ui.column(
                    4,
                    ui.input_numeric('rolling_avg', 'Rolling average', value=20, min=1)
                ),
                ui.column(
                    4,
                    ui.input_selectize('stat_color', 'Line color', ['indigo', 'darkblue', 'grey', 'black', 'burlywood'], selected='grey')
                )
            ),
            ui.row(
                ui.h6('Second plot:', class_='section-title'),
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
                ui.h6('Third plot:',  class_='section-title'),
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
                ),
                ui.column(
                    4,
                ui.h6('Annotation plot'),
                    ui.input_selectize('feature_display', 'Feature to display', ['None']),
                    ui.input_selectize('feature_color', 'Feature color', ['indigo', 'darkblue', 'grey', 'black', 'burlywood'], selected='grey')
                ),
            )
        ),
        ui.nav_panel(
            'Visualization',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_switch('seq_names', 'show sequence names', value=False),
                    ui.input_radio_buttons('stat_type', ui.h6('Select plot type:'), ['gc', 'entropy', 'coverage', 'identity'], selected='gc'),
                    ui.input_radio_buttons( 'alignment_type', ui.h6('Select plot type:'), ['identity', 'similarity'], selected='identity'),
                    ui.input_radio_buttons('annotation', ui.h6('Select plot type:'), ['SNPs','Conserved ORFs', 'Annotation'], selected='Annotation')
                ),
                ui.output_plot('msa_plot', height='90vh', width='90vw'),
                ),
                ui.div(
                    ui.input_slider('zoom_range', 'Zoom Range (Start - Stop):', min=0, max=1000, value=(0, 1000), step=1),
                    style='position: absolute; top: 5px; right: 10px;'
            )
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
            ui.update_slider('zoom_range', max=alignment_length, value=(0, alignment_length))

    @reactive.Effect
    @reactive.event(input.annotation_file)
    def load_annotation():
        annotation_file = input.annotation_file()
        if annotation_file:
            ann = explore.Annotation(reactive.alignment.get(), annotation_file[0]['datapath'])
            reactive.annotation.set(ann)

    @reactive.Effect
    @reactive.event(input.alignment_file)
    def update_reference():
        aln = reactive.alignment.get()
        ui.update_selectize(
            id='reference', choices=['first' ,'consensus']+list(aln.alignment.keys()), selected='first'
        )

    @reactive.Effect
    @reactive.event(input.alignment_file)
    def subsitution_matrix():
        aln = reactive.alignment.get()
        ui.update_selectize(
            id='matrix',
            choices=list(config.SUBS_MATRICES[aln.aln_type].keys()),
            selected='BLOSUM65' if aln.aln_type == 'AA' else 'TRANS',
        )

    @reactive.Effect
    @reactive.event(input.annotation_file)
    def define_annotation_to_display():
        ann = reactive.annotation.get()
        ui.update_selectize(
            id='feature_display',
            choices=list(ann.features.keys()),
            selected=list(ann.features.keys())[0]
        )

    @output
    @render.plot
    def msa_plot():
        # Ensure an alignment file has been uploaded
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()
        if not aln:
            return None

        # set the reference sequence
        if input.reference() == 'first':
            aln.reference_id = list(aln.alignment.keys())[0]
        elif input.reference() == 'consensus':
            aln.reference_id = None
        else:
            aln.reference_id = input.reference()

        # Update zoom level from slider
        aln.zoom = tuple(input.zoom_range())

        # Adjust the size depending on the number of alignment sequences
        aln_len, current_height_ratio = len(aln.alignment.keys()), 0.1
        for n_seq, height_ratio in zip(
                [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500],
                [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 7.5, 10]
        ):
            if aln_len >= n_seq:
                current_height_ratio = height_ratio
                break
        # Prepare the plot with 3 subplots
        fig, axes = plt.subplots(nrows=3, height_ratios=[0.05, current_height_ratio, 0.1])

        # Subplot 1: Stats Plot
        draw.stat_plot(
            aln,
            axes[0],
            stat_type=input.stat_type(),
            line_width=1,
            rolling_average=input.rolling_avg(),
            line_color=input.stat_color()
        )

        # Subplot 2: Alignment Plot (Identity or Similarity)
        if input.alignment_type() == 'identity':
            draw.identity_alignment(
                aln, axes[1],
                fancy_gaps=input.show_gaps(),
                show_gaps=input.show_gaps(),
                show_mask=True,
                show_mismatches=True,
                show_ambiguities=True,
                reference_color='lightsteelblue',
                show_seq_names=input.seq_names(),
                show_x_label=True if input.annotation() == 'Annotation' and not input.annotation_file() else False,
                show_legend=input.show_legend()
            )
        else:
            draw.similarity_alignment(
                aln, axes[1],
                fancy_gaps=input.show_gaps(),
                show_gaps=input.show_gaps(),
                matrix_type=input.matrix(),
                show_seq_names=input.seq_names(),
                show_cbar=input.show_legend(),
                cbar_fraction=0.02,
                show_x_label=True if input.annotation() == 'Annotation' and not input.annotation_file() else False
            )

        # Subplot 3: Annotation or ORF Plot
        if input.annotation() == 'Annotation' and input.annotation_file():
            draw.annotation_plot(
                aln, ann,
                axes[2],
                feature_to_plot=input.feature_display(),
                show_x_label=True
            )
        elif input.annotation() == 'Conserved ORFs':
            draw.orf_plot(
                aln, axes[2],
                cmap='hsv',
                non_overlapping_orfs=input.non_overlapping(),
                show_x_label=True,
                show_cbar=True,
                cbar_fraction=0.2,
                min_length=input.min_orf_length()
            )
        elif input.annotation() == 'SNPs':
            draw.variant_plot(
                aln, axes[2],
                show_x_label=True,
                lollisize=(input.stem_size(), input.head_size()),
                show_legend=input.show_legend_variants()
            )
        # turn everything off
        else:
            axes[2].axis('off')

        return fig


app = App(app_ui, server)
