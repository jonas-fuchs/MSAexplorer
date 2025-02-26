from shiny import App, render, ui, reactive
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

# Define the UI
app_ui = ui.page_fluid(
    ui.h2('MSAexplorer'),
    ui.navset_tab(
        ui.nav_panel(
            'Upload Files',
            ui.input_file('alignment_file', 'Upload Alignment File (.aln)', multiple=False),
            ui.input_file('annotation_file', 'Optional Annotation File (.gff3, .bed, .gb)', multiple=False),
        ),
        ui.nav_panel(
            'Advanced Settings',
            ui.input_numeric('rolling_avg', 'Rolling Average (Statistic Plot):', value=20, min=1),
            ui.input_checkbox('show_gaps', 'Show Gaps in Alignment:', value=True),
            ui.input_checkbox('show_annotation', 'Show Annotations (else ORF Plot):', value=True),
            ui.input_numeric('min_orf_length', 'Minimum ORF Length:', value=150, min=1),
            ui.input_selectize('reference', 'Reference sequence', ['first' ,'consensus'], selected='first'),
        ),
        ui.nav_panel(
            'Visualization',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_radio_buttons('stat_type', 'Select Statistic Type:', ['gc', 'entropy', 'coverage', 'identity'], selected='gc'),
                    ui.input_radio_buttons( 'alignment_type', 'Alignment Type:', ['identity', 'similarity'], selected='identity'),
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
    @reactive.event(input.alignment_file)
    def update_reference():
        aln = reactive.alignment.get()
        ui.update_selectize(
            id='reference', choices=['first' ,'consensus']+list(aln.alignment.keys()), selected='first'
        )

    @output
    @render.plot
    def msa_plot():
        # Ensure an alignment file has been uploaded
        aln = reactive.alignment.get()
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
            rolling_average=input.rolling_avg(),
            line_color='black'
        )

        # Subplot 2: Alignment Plot (Identity or Similarity)
        if input.alignment_type() == 'identity':
            draw.identity_alignment(
                aln, axes[1],
                show_gaps=input.show_gaps(),
                show_mask=True,
                show_mismatches=True,
                reference_color='lightsteelblue',
                show_seq_names=False,
                show_x_label=False,
                show_legend=True
            )
        else:
            draw.similarity_alignment(
                aln, axes[1],
                fancy_gaps=True,
                show_gaps=input.show_gaps(),
                matrix_type='TRANS',
                show_cbar=True,
                cbar_fraction=0.02,
                show_x_label=False
            )

        # Subplot 3: Annotation or ORF Plot
        if input.show_annotation() and input.annotation_file():
            annotation_file = input.annotation_file()
            draw.annotation_plot(
                aln, annotation_file[0]['datapath'],
                axes[2],
                feature_to_plot='gene',
                show_x_label=True
            )
        else:
            draw.orf_plot(
                aln, axes[2],
                cmap='hsv',
                non_overlapping_orfs=False,
                show_cbar=True,
                cbar_fraction=0.2,
                min_length=input.min_orf_length()
            )

        # Adjust layout
        return fig


app = App(app_ui, server)
