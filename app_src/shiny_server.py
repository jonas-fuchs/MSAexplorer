"""
This module creates the Server logic
"""

# built-in
import tempfile

# libs
from shiny import render, ui, reactive
import matplotlib.pyplot as plt
from matplotlib import colormaps

# load in app resources
from app_src.shiny_plots import set_aln, create_msa_plot, create_analysis_custom_heatmap, create_freq_heatmap
from shinywidgets import render_widget

# msaexplorer
from msaexplorer import explore, config, export


def server(input, output, session):
    """
    creates the server logic
    """

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
    def prepare_minimal_inputs(zoom: bool = True, window_size: bool = False, ref: bool = False, left_plot: bool = False, right_plot: bool = False):
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
        # left analysis plot
        if left_plot:
            inputs['analysis_plot_type_left'] = input.analysis_plot_type_left()
            if inputs['analysis_plot_type_left'] != 'Off':
                inputs['additional_analysis_options_left'] = input.additional_analysis_options_left()
        # right analysis plot
        if right_plot:
            inputs['analysis_plot_type_right'] = input.analysis_plot_type_right()
            if inputs['analysis_plot_type_right'] != 'Off':
                inputs['additional_analysis_options_right'] = input.additional_analysis_options_right()
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
            ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Annotation'], selected='Off')
        else:
            ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Conserved ORFs', 'Annotation'],
                                selected='Off')

    # Updates the plot container
    @reactive.Effect
    async def update_height():
        """
        update the plot container height -> Sends a message that is picked up by the js and updates the CSS height
        property for the plot container
        """
        new_plot_height = f'{input.increase_height() * 100}vh'

        await session.send_custom_message("update-plot-container-height", {'height': new_plot_height})

    #TODO: Also try to conditionally exclude the Annotation column

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

                # add specific settings to settings tab
                # on each upload first remove potential ui elements
                ui.remove_ui(selector="#orf_column")
                # then add the ui element
                if aln.aln_type != 'AA':
                    ui.insert_ui(
                        ui.column(
                            4,
                            ui.h6('ORF plot'),
                            ui.input_numeric('min_orf_length', 'Length', value=150, min=1),
                            ui.input_selectize('color_mapping', 'Colormap ORF identity', choices=list(colormaps.keys()),
                                               selected='jet'),
                            ui.input_switch('non_overlapping', 'non-overlapping', value=False),
                            id='orf_column'
                        ),
                        selector='#snp_column',
                        where='afterEnd'
                    )

                # Update reference
                for id in ['reference', 'reference_2']:
                    ui.update_selectize(
                        id=id, choices=['first', 'consensus'] + list(aln.alignment.keys()), selected='first'
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
        load the annotation - catches errors if its the wrong format and displays it (otherwise app_src would crash)
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

        # remove prior ui elements and insert specific once prior the download format div
        ui.remove_ui(selector="div:has(> #download_type_options_1-label)")
        ui.remove_ui(selector="div:has(> #download_type_options_2-label)")
        ui.remove_ui(selector="div:has(> #reference_2-label)")
        if input.download_type() == 'SNPs':
            ui.update_selectize('download_format', choices=['vcf', 'tabular'])
            ui.insert_ui(
                ui.input_selectize('download_type_options_1', label='include ambiguous snps', choices=['Yes', 'No'], selected='No'),
                selector='#download_format-label',
                where='beforeBegin'
            )
            ui.insert_ui(
                ui.input_selectize('reference_2', 'Reference', ['first', 'consensus'], selected='first'),
                selector='#download_format-label',
                where='beforeBegin'
            )
        elif input.download_type() == 'consensus':
            ui.update_selectize('download_format', choices=['fasta'])
            ui.insert_ui(
                ui.input_selectize('download_type_options_1', label='Use ambiguous characters (only nt)', choices=['Yes', 'No'], selected='No'),
                selector='#download_format-label',
                where='beforeBegin'
            )
            ui.insert_ui(
                ui.input_numeric('download_type_options_2', label='Frequency threshold', value=0, min=0, max=1, step=0.1),
                selector='#download_format-label',
                where='beforeBegin'
            )

    # TODO: Download distance matrices
    # TODO: Download ORFs as bed
    # TODO: Download entropy
    # TODO: Download gc
    # TODO: Download ts/tv
    # TODO: Download coverage
    # TODO: Download similarity
    # TODO: Download identity
    # TODO: Download char frequencies
    # TODO: Download reverse complement
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
                if input.reference_2() == 'first':
                    aln.reference_id = list(aln.alignment.keys())[0]
                elif input.reference_2()  == 'consensus':
                    aln.reference_id = None
                else:
                    aln.reference_id = input.reference_2()
                download_data = export.snps(aln.get_snps(include_ambig=True if input.download_type_options_1 == 'Yes' else False), format_type=download_format)
                prefix = 'SNPs_'
            elif input.download_type() == 'consensus':
                if input.download_type_options_1() == 'Yes' and aln.aln_type != 'AA':
                    download_data = export.fasta(aln.get_consensus(threshold=input.download_type_options_2(), use_ambig_nt=True))
                else:
                    download_data = export.fasta(aln.get_consensus(threshold=input.download_type_options_2()))
                prefix = 'consensus_'

            # Create a temporary file for the download
            with tempfile.NamedTemporaryFile(prefix=prefix, suffix=f'.{download_format}', delete=False) as tmpfile:
                tmpfile.write(download_data.encode('utf-8'))
                tmpfile.flush()  # Ensure data is written to disk

                return tmpfile.name

        # send a notification to the user if no alignment was uploaded
        except FileNotFoundError:
            ui.notification_show(ui.tags.div(
                'No alignment was uploaded.',
                style="color: red; font-weight: bold;"
            ), duration=10)

            return None
        # or if the Threshold was not set correctly
        except ValueError:
            ui.notification_show(ui.tags.div(
                'Threshold frequency value has to be between 0 and 1.',
                style="color: red; font-weight: bold;"
            ), duration=10)

            return None

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

    #TODO: Conditionally exclude options

    @reactive.Effect
    @reactive.event(input.analysis_plot_type_left)
    def update_additional_options_left():
        """
        Update UI for the left plot in the analysis tab
        """
        # ensure that it is switched back
        if input.analysis_plot_type_left() == 'Off':
            ui.update_selectize('additional_analysis_options_left', choices=['None'], selected='None')
        if input.analysis_plot_type_left() == 'Pairwise identity':
            ui.update_selectize(
                'additional_analysis_options_left',
                choices={
                    'ghd': 'global hamming distance',
                    'lhd': 'local hamming distance',
                    'ged': 'gap excluded distance',
                    'gcd': 'gap compressed distance'
                },
                selected='ghd'
            )

    @reactive.Effect
    @reactive.event(input.analysis_plot_type_right)
    def update_additional_options_right():
        """
        Update UI for the left plot in the analysis tab
        """
        # ensure that it is switched back
        if input.analysis_plot_type_right() == 'Off':
            ui.update_selectize('additional_analysis_options_right', choices=['None'], selected='None')

    @render.text
    def analysis_info_left():
        """
        show custom text for additional options
        """
        selected_option = input.additional_analysis_options_left()
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

    @render.text
    def analysis_info_right():
        """
        show custom text for additional options
        """
        return None


    @render_widget
    def analysis_custom_heatmap():
        """
        Create the heatmap
        """
        aln = reactive.alignment.get()

        return create_analysis_custom_heatmap(aln, prepare_minimal_inputs(left_plot=True, window_size=True))


    @render_widget
    def analysis_char_freq_heatmap():
        """
        Create character frequency heatmap
        """
        aln = reactive.alignment.get()

        return create_freq_heatmap(aln,  prepare_minimal_inputs(right_plot=True, window_size=True))

