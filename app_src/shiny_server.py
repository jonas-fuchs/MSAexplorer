"""
This module creates the Server logic
"""

# built-in
import tempfile
from typing import Callable, Dict
import asyncio
import urllib.request

# libs
import numpy as np
from numpy import ndarray
from shiny import render, ui, reactive
import matplotlib.pyplot as plt
from matplotlib import colormaps

# try to load pyfamsa and pytrimal for static website compatibility
try:
    from pyfamsa import Aligner, Sequence
    pyfamsa_check = True
except ImportError:
    pyfamsa_check = False

try:
    from pytrimal import Alignment, AutomaticTrimmer
    pytrimal_check = True
except ImportError:
    pytrimal_check = False

# load in app resources
from app_src.shiny_plots import set_aln, create_msa_plot, create_analysis_custom_heatmap, create_freq_heatmap, create_recovery_heatmap
from shinywidgets import render_widget

# msaexplorer
from msaexplorer import explore, config, export, draw


def server(input, output, session):
    """
    creates the server logic
    """

    # ini several params
    reactive.alignment = reactive.Value(None)
    reactive.annotation = reactive.Value(None)
    updating_from_slider = False
    updating_from_numeric = False

    # remove either pytrimal button or pyfamsa button if not installed
    if not pytrimal_check:
        ui.remove_ui(selector="#trim")
    elif not pyfamsa_check:
        ui.remove_ui(selector="#align")


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
        inputs['example_alignment'] = input.example_alignment()

        # STATISTICS (first plot)
        if inputs['stat_type'] not in ['Off', 'sequence logo']:
            inputs['plot_1_size'] = input.plot_1_size()
            inputs['rolling_average'] = input.rolling_avg()
            inputs['stat_color'] = input.stat_color()
        elif inputs['stat_type'] == 'sequence logo':
            inputs['plot_1_size'] = input.plot_1_size()
            inputs['logo_coloring'] = input.logo_coloring()
            relative_msa_width = window_width / (inputs['zoom_range'][1]-inputs['zoom_range'][0])
            if relative_msa_width >= 11:
                inputs['logo_type'] = 'logo'
            else:
                inputs['logo_type'] = 'stacked'

        # ALIGNMENT (second plot)
        if inputs['alignment_type'] != 'Off':
            inputs['reference_color'] = input.reference_color()
            inputs['plot_2_size'] = input.plot_2_size()
            inputs['fancy_gaps'] = input.fancy_gaps()
            inputs['show_mask'] = input.show_mask()
            inputs['show_gaps'] = input.show_gaps()
            inputs['show_legend'] = input.show_legend()
            inputs['show_ambiguities'] = input.show_ambiguities()
            inputs['identity_coloring'] = input.identity_coloring()
            if inputs['alignment_type'] == 'similarity':
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
                inputs['show_sequence_all'] = input.show_sequence_all()
            else:
                inputs['show_sequence'] = False
                inputs['show_sequence_all'] = False
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
                inputs['snp_coloring'] = input.snp_coloring()
                inputs['stem_size'] = input.stem_size()
                inputs['head_size'] = input.head_size()
                inputs['show_legend_third_plot'] = input.show_legend_third_plot()

        return inputs


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
        if window_size:
            inputs['dimensions'] = input.window_dimensions()

        return inputs


    def read_in_annotation(annotation_file:str):
        """
        Read in an annotation and update the ui accordingly
        """
        try:
            ann = explore.Annotation(reactive.alignment.get(), annotation_file[0]['datapath'])
        # when loading example data there is no datapath
        except TypeError:
            ann = explore.Annotation(reactive.alignment.get(), annotation_file)
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


    def is_sequence_list(filepath:str):
        """
        Quickly check if the file is a list of sequences, not a full file
        """
        with open(filepath) as f:
            lines = [line.strip() for line in f if not line.startswith(">")]
        lengths = set(len(line) for line in lines if line)
        return len(lengths) > 1


    def align_sequences(alignment_file, n_threads:int, guide_tree:str, refine:bool, keep_duplicates:bool):
        """
        What happens when the align sequences button is pressed.
        """
        sequences, seq_id = [], None
        with open(alignment_file[0]['datapath'], 'r') as file:
            for i, line in enumerate(file):
                line = line.strip()
                if line.startswith(">"):
                    if seq_id is not None:
                        seq = ''.join(seq).upper()
                        sequences.append(Sequence(seq_id.encode(), seq.encode()))  # the wrapper needs bit encoding
                    seq_id = line[1:]
                    seq = []
                else:
                    seq.append(line)
            # append the last sequence
            if seq_id is not None and seq:
                seq = ''.join(seq).upper()
                sequences.append(Sequence(seq_id.encode(), seq.encode()))
            # define the aligner method
            aligner = Aligner(threads=n_threads, guide_tree=guide_tree, refine=refine, keep_duplicates=keep_duplicates)
            # align sequences
            alignment_new = aligner.align(sequences)
            aln_string = ''
            for sequence in alignment_new:
                aln_string = f"{aln_string}>{sequence.id.decode()}\n{sequence.sequence.decode()}\n"  # create str and decode

        return aln_string


    def transfer_alignment_to_trimAI(aln):
        """Compatibility between trimAI and msaexplorer by reading and writing to a temporary file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tmp', delete=False) as tmp:
            tmp.write(export.fasta(aln._alignment))
            tmp.flush()
            return Alignment.load(tmp.name)


    def trim_sequences(aln, method:str):
        """
        What happens when the trim button is pressed.
        """

        aln = transfer_alignment_to_trimAI(aln)
        trimmer = AutomaticTrimmer(method=method)
        trimmed_aln = trimmer.trim(aln)
        aln_string = ''
        for name, seq in zip(trimmed_aln.names, trimmed_aln.sequences):
            aln_string = f"{aln_string}>{name.decode()}\n{seq}\n"

        return aln_string


    def finalize_loaded_alignment(aln, annotation_file):
        """
        All the necessary steps to load an alignment into the app are done here.
        """

        aln.reference_id = list(aln.alignment.keys())[0]
        reactive.alignment.set(aln)

        alignment_length = len(next(iter(aln.alignment.values()))) - 1
        ui.update_slider('zoom_range', max=alignment_length - 1, value=(0, alignment_length - 1))

        ui.remove_ui(selector="#orf_column")
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

        for id in ['reference', 'reference_2']:
            ui.update_selectize(id=id, choices=['first', 'consensus'] + list(aln.alignment.keys()), selected='first')

        ui.update_selectize(
            id='matrix',
            choices=list(config.SUBS_MATRICES[aln.aln_type].keys()),
            selected='BLOSUM65' if aln.aln_type == 'AA' else 'TRANS',
        )

        aln_len, seq_threshold = len(aln.alignment.keys()), 5
        for ratio in config.STANDARD_HEIGHT_RATIOS.keys():
            if aln_len >= ratio:
                seq_threshold = ratio

        ui.update_numeric('plot_1_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][0])
        ui.update_numeric('plot_2_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][1])
        ui.update_numeric('plot_3_size', value=config.STANDARD_HEIGHT_RATIOS[seq_threshold][2])

        if aln.aln_type == 'AA':
            ui.update_selectize('stat_type',
                                choices=['Off', 'sequence logo', 'entropy', 'coverage', 'identity', 'similarity'],
                                selected='Off')
            ui.update_selectize('download_type',
                                choices=['alignment','SNPs', 'consensus', 'character frequencies', '% recovery', 'entropy',
                                         'coverage', 'mean identity', 'mean similarity'])
            ui.update_selectize('annotation', choices=['Off', 'SNPs'])
            ui.update_selectize('identity_coloring', choices=['None', 'standard', 'clustal', 'zappo', 'hydrophobicity'])
            ui.update_selectize('logo_coloring', choices=['standard', 'clustal', 'zappo', 'hydrophobicity'])
            ui.update_selectize('snp_coloring', choices=['standard', 'clustal', 'zappo', 'hydrophobicity'])
        else:
            ui.update_selectize('stat_type',
                                choices=['Off', 'sequence logo', 'gc', 'entropy', 'coverage', 'identity', 'similarity',
                                         'ts tv score'], selected='Off')
            ui.update_selectize('download_type', choices=['alignment', 'SNPs', 'consensus', 'character frequencies', '% recovery',
                                                          'reverse complement alignment', 'conserved orfs', 'gc',
                                                          'entropy', 'coverage', 'mean identity', 'mean similarity',
                                                          'ts tv score'])
            ui.update_selectize('annotation', choices=['Off', 'SNPs', 'Conserved ORFs'])
            ui.update_selectize('identity_coloring',
                                choices=['None', 'standard', 'standard', 'purine_pyrimidine', 'strong_weak'])
            ui.update_selectize('logo_coloring', choices=['standard', 'standard', 'purine_pyrimidine', 'strong_weak'])
            ui.update_selectize('snp_coloring', choices=['standard', 'standard', 'purine_pyrimidine', 'strong_weak'])

        if annotation_file:
            read_in_annotation(annotation_file)


    def show_alignment_error(e, alignment_file):
        """
        Triggers the ui notification to show an error message.
        """
        ui.notification_show(ui.tags.div(f'Error: {e}', style="color: red; font-weight: bold;"), duration=10)
        if alignment_file and is_sequence_list(alignment_file[0]['datapath']) and pyfamsa_check:
            ui.notification_show(ui.tags.div(
                'It seems like you have uploaded a list of sequences. No worries, go ahead and align them in the app.',
                style="color: black"
            ), duration=10)


    ##### handel everything upload and sequence processing related #####
    @ui.bind_task_button(button_id='trim')
    @reactive.extended_task
    async def trimming(aln, method):
        """
        Asynchronous task handler for trimming sequences. This function binds a button
        action to a reactive extended task for sequence trimming.
        """
        return await asyncio.to_thread(trim_sequences, aln, method)

    @ui.bind_task_button(button_id='align')
    @reactive.extended_task
    async def aligning(alignment_file, n_threads:int, guide_tree:str, refine:bool, keep_duplicates:bool):
        """
        Asynchronous task handler for aligning sequences. This function binds a button
        action to a reactive extended task for sequence alignment.
        """
        return await asyncio.to_thread(align_sequences, alignment_file, n_threads, guide_tree, refine, keep_duplicates)

    @reactive.Effect
    @reactive.event(input.align)
    def alignment_task():
        """
        Handles the alignment execution.
        """
        alignment_file = input.alignment_file()
        if alignment_file:
            # decode guide tree options
            guide_tree_map = {
                'MST + Prim single linkage': 'sl',
                'SLINK single linkage': 'slink',
                'UPGMA': 'upgma',
                'neighbour joining': 'nj'
            }
            aligning.invoke(
                alignment_file=alignment_file,
                n_threads=input.n_threads(),
                guide_tree=guide_tree_map[input.guide_tree()],
                refine=input.refine(),
                keep_duplicates=input.keep_duplicates()
            )
            ui.notification_show(ui.tags.div(
                'Alignment execution started.',
                style="color: black; font-weight: bold;"
            ), duration=10)

    @reactive.Effect
    @reactive.event(input.trim)
    def trim_task():
        """
        Handles the trimming execution.
        """
        aln = reactive.alignment.get()
        if aln is not None:
            trimming.invoke(
                aln=aln,
                method=input.trim_method()
            )
            ui.notification_show(ui.tags.div(
                'Trimming execution started.',
                style="color: black; font-weight: bold;"
            ), duration=10)

    @reactive.effect
    @reactive.event(input.cancel)
    def execution_cancel():
        """
        Handles the cancel button.
        """
        aligning.cancel()
        trimming.cancel()
        ui.notification_show(ui.tags.div(
            'Execution cancelled.',
            style="color: black; font-weight: bold;"
        ), duration=10)

    @reactive.Effect
    @reactive.event(aligning.result)
    def load_aligned_result():
        """
        The function listens for the completion of the `align_sequences_task` and loads the aligned
        result into an MSA (Multiple Sequence Alignment) object. This function is reactive, meaning
        it automatically responds to changes in the `align_sequences_task.result` value. It also
        handles any potential errors during the alignment loading process.
        """
        if not pyfamsa_check:
            return None
        alignment_file = input.alignment_file()
        annotation_file = input.annotation_file()
        try:
            alignment_finished = aligning.result()
            ui.notification_show(ui.tags.div(
                'Alignment was successful.',
                style="color: green; font-weight: bold;"
            ), duration=10)
            aln = explore.MSA(alignment_finished, reference_id=None, zoom_range=None)
            finalize_loaded_alignment(aln, annotation_file)
        except Exception as e:
            show_alignment_error(e, alignment_file)

    @reactive.Effect
    @reactive.event(trimming.result)
    def load_trimmed_result():
        """
        reload trimmed alignment
        """
        if not pytrimal_check:
            return None
        annotation_file = input.annotation_file()
        try:
            alignment_finished = trimming.result()
            ui.notification_show(ui.tags.div(
                'Trimming was successful.',
                style="color: green; font-weight: bold;"
            ), duration=10)
            aln = explore.MSA(alignment_finished, reference_id=None, zoom_range=None)
            finalize_loaded_alignment(aln, annotation_file)
        except Exception as e:
            ui.notification_show(ui.tags.div(
                f'{e}',
                style="color: red"
            ), duration=10)

    @reactive.Effect
    @reactive.event(input.alignment_file)
    def load_direct_alignment():
        """
        Monitors if an new alignment file is uploaded and loads it directly into the app.
        """
        alignment_file = input.alignment_file()
        annotation_file = input.annotation_file()
        try:
            aln = explore.MSA(alignment_file[0]['datapath'], reference_id=None, zoom_range=None)
            finalize_loaded_alignment(aln, annotation_file)
        except Exception as e:
            show_alignment_error(e, alignment_file)

    @reactive.Effect
    @reactive.event(input.example_alignment)
    def load_example(alignment_file=True):
        try:
            with urllib.request.urlopen('https://raw.githubusercontent.com/jonas-fuchs/MSAexplorer/refs/heads/master/example_alignments/DNA.fasta') as resp:
                alignment_file = resp.read().decode("utf-8")
            with urllib.request.urlopen('https://raw.githubusercontent.com/jonas-fuchs/MSAexplorer/refs/heads/master/example_alignments/DNA_RNA.gb') as resp:
                annotation_file = resp.read().decode("utf-8")
            aln = explore.MSA(alignment_file)
            finalize_loaded_alignment(aln, annotation_file)
            ui.notification_show("Example alignment and annotation file loaded.")
        except ValueError as e:
            print(e)
            ui.notification_show(ui.tags.div(
                f'Error: {e}',
                style="color: red; font-weight: bold;"
            ), duration=10)

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


    #### Handle everything plotting related ####
    @reactive.Effect
    async def update_height():
        """
        update the plot container height -> Sends a message that is picked up by the js and updates the CSS height
        property for the plot container
        """
        new_plot_height = f'{input.increase_height() * 100}vh'

        await session.send_custom_message("update-plot-container-height", {'height': new_plot_height})

    @reactive.effect
    def update_zoom_boxes():
        """
        Update the zoom boxes when the slider changes.
        """
        nonlocal updating_from_slider
        # do not update if we are already updating from a slider change or from a numeric change.
        # importantly, this now includes updating from slider check which prohibits a too early update
        # of when trying to update the zoom again via the slider before the previous plot was generated.
        if updating_from_numeric or updating_from_slider:
            return
        updating_from_slider = True
        zoom_range = input.zoom_range()
        ui.update_numeric("zoom_start", value=zoom_range[0])
        ui.update_numeric("zoom_end", value=zoom_range[1])
        updating_from_slider = False

    @reactive.effect
    def update_zoom_slider():
        """
        Update the zoom slider when the zoom input changes.
        """
        nonlocal updating_from_numeric
        if updating_from_numeric or updating_from_slider:
            return
        updating_from_numeric = True
        start, end = input.zoom_start(), input.zoom_end()
        # make sure set values make sense
        if start is None or end is None:
            return
        if start >= end:
            end = start +1
        ui.update_slider("zoom_range", value=(start, end))
        updating_from_numeric = False

    @output
    @render.plot
    def msa_plot():
        """
        plot the alignment
        """
        aln = reactive.alignment.get()
        ann = reactive.annotation.get()

        return create_msa_plot(aln, ann, prepare_inputs())


    #### Handle everything download related ###
    @reactive.Effect
    @reactive.event(input.download_type)
    def update_download_options():
        """
        Update the UI - which download options are displayed
        """

        # remove prior ui elements and insert specific once prior the download format div
        aln = reactive.alignment.get()
        ui.remove_ui(selector="div:has(> #download_type_options_1-label)")
        ui.remove_ui(selector="div:has(> #download_type_options_2-label)")
        ui.remove_ui(selector="div:has(> #download_type_options_3-label)")
        ui.remove_ui(selector="div:has(> #reference_2-label)")
        if input.download_type() == 'SNPs' or input.download_type() == '% recovery':
            if input.download_type() == '% recovery':
                ui.update_selectize('download_format', choices=['csv', 'tabular'])
            else:
                ui.update_selectize('download_format', choices=['vcf', 'tabular'])
            if input.download_type() == 'SNPs':
                ui.insert_ui(
                    ui.input_selectize('download_type_options_1', label='include ambiguous snps', choices=['Yes', 'No'], selected='No'),
                    selector='#download_format-label',
                    where='beforeBegin'
                )
            ui.insert_ui(
                ui.input_selectize(
                    'reference_2', 'Reference', ['first', 'consensus'], selected='first'
                ),
                selector='#download_format-label',
                where='beforeBegin'
            ) if aln is None else ui.insert_ui(
                ui.input_selectize(
                    id='reference_2', label='Reference', choices=['first', 'consensus'] + list(aln.alignment.keys()), selected='first'
                ),
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
        elif input.download_type() in ['entropy', 'coverage', 'mean identity', 'mean similarity', 'ts tv score', 'gc']:
            ui.update_selectize('download_format', choices=['csv', 'tabular'])
            ui.insert_ui(
                ui.input_numeric('download_type_options_1', label='Rolling average', value=1, min=1, step=1),
                selector='#download_format-label',
                where='beforeBegin'
            )
        elif input.download_type() == 'conserved orfs':
            ui.update_selectize('download_format', choices=['bed'])
            ui.insert_ui(
                ui.input_numeric('download_type_options_1', label='min length', value=100, min=6, step=1),
                selector='#download_format-label',
                where='beforeBegin'
            )
            ui.insert_ui(
                ui.input_numeric('download_type_options_2', label='identity cutoff', value=95, min=0, max=100, step=1),
                selector='#download_format-label',
                where='beforeBegin'
            )
            ui.insert_ui(
                ui.input_selectize('download_type_options_3', label='Allow overlapping orfs?', choices=['Yes', 'No'], selected='Yes'),
                selector='#download_format-label',
                where='beforeBegin'
            )
        elif input.download_type() in ['alignment', 'reverse complement alignment']:
            ui.update_selectize('download_format', choices=['fasta'])
        elif input.download_type() == 'character frequencies':
            ui.update_selectize('download_format', choices=['tabular', 'csv'])

    @render.download()
    def download_stats():
        """
        Download various files in standard format
        """

        # helper functions
        def _snp_option():
            if input.reference_2() == 'first':
                aln.reference_id = list(aln.alignment.keys())[0]
            elif input.reference_2() == 'consensus':
                aln.reference_id = None
            else:
                aln.reference_id = input.reference_2()
            download_data = export.snps(
                aln.get_snps(include_ambig=True if input.download_type_options_1 == 'Yes' else False),
                format_type=download_format)

            return download_data, 'SNPs_'

        def _consensus_option():
            if input.download_type_options_1() == 'Yes' and aln.aln_type != 'AA':
                download_data = export.fasta(
                    sequence=aln.get_consensus(threshold=input.download_type_options_2(), use_ambig_nt=True),
                    header='ambiguous_consensus',
                )
            else:
                download_data = export.fasta(
                    sequence=aln.get_consensus(threshold=input.download_type_options_2()),
                    header='consensus',
                )

            return download_data, 'consensus_'

        def _stat_option():
            # create function mapping
            stat_functions: Dict[str, Callable[[], list | ndarray]] = {
                'gc': aln.calc_gc,
                'entropy': aln.calc_entropy,
                'coverage': aln.calc_coverage,
                'mean identity': aln.calc_identity_alignment,
                'mean similarity': aln.calc_similarity_alignment,
                'ts tv score': aln.calc_transition_transversion_score
            }
            # raise error for rolling average
            if input.download_type_options_1() < 1 or input.download_type_options_1() > aln.length:
                raise ValueError('Rolling_average must be between 1 and length of alignment.')
            # define seperator
            seperator = '\t' if input.download_format() == 'tabular' else ','
            # define which stat type to exprt
            for stat_type in ['entropy', 'mean similarity', 'coverage', 'mean identity', 'ts tv score', 'gc']:
                if stat_type == input.download_type():
                    break
            # use correct function
            data = stat_functions[stat_type]()
            # calculate the mean (identical to draw module of msaexplorer)
            if stat_type in ['mean identity', 'mean similarity']:
                # for the mean nan values get handled as the lowest possible number in the matrix
                data = np.nan_to_num(data, True, -1 if stat_type == 'identity' else 0)
                data = np.mean(data, axis=0)
            # apply rolling average
            data = draw._moving_average(data, input.download_type_options_1(), None, aln.length)[0]
            # create download data
            download_data = export.stats(stat_data=data, seperator=seperator)
            prefix = f'{stat_type}_'.replace(' ', '_')  # sanitize prefix

            return download_data, prefix

        def _orf_option():
            if input.download_type_options_3() == 'Yes':
                data = aln.get_conserved_orfs(min_length=input.download_type_options_1(), identity_cutoff=input.download_type_options_2())
            else:
                data = aln.get_non_overlapping_conserved_orfs(min_length=input.download_type_options_1(),identity_cutoff=input.download_type_options_2())

            return export.orf(data, aln.reference_id.split(' ')[0]), 'orfs_' if input.download_type_options_3() == 'Yes' else 'non_overlapping_conserved_orfs_'

        def _alignment_option():
            data = aln.alignment

            return export.fasta(sequence=data), 'alignment_'

        def _reverse_complement_option():
            data = aln.calc_reverse_complement_alignment()

            return export.fasta(sequence=data), 'rc_alignment_'

        def _char_freq_option():
            data = aln.calc_character_frequencies()

            return export.character_freq(data, seperator='\t' if input.download_format() == 'tabular' else ','), 'char_freq_'

        def _percent_recovery_option():
            if input.reference_2() == 'first':
                aln.reference_id = list(aln.alignment.keys())[0]
            elif input.reference_2() == 'consensus':
                aln.reference_id = None
            else:
                aln.reference_id = input.reference_2()

            download_data = export.percent_recovery(
                aln.calc_percent_recovery(),
                seperator='\t' if input.download_format() == 'tabular' else ',')

            return download_data, 'recovery_'

        try:
            # Initialize
            download_format = input.download_format()
            aln = reactive.alignment.get()
            if aln is None:
                raise FileNotFoundError("No alignment data available. Please upload an alignment.")
            # Create download data for SNPs
            if input.download_type() == 'SNPs':
                export_data = _snp_option()
            # Create download data for consensus
            elif input.download_type() == 'consensus':
                export_data = _consensus_option()
            # Create download data for stats
            elif input.download_type() in ['entropy', 'mean similarity', 'coverage', 'mean identity', 'ts tv score', 'gc']:
                export_data = _stat_option()
            elif input.download_type() == 'conserved orfs':
                export_data = _orf_option()
            elif input.download_type() == 'reverse complement alignment':
                export_data = _reverse_complement_option()
            elif input.download_type() == 'alignment':
                export_data = _alignment_option()
            elif input.download_type() == 'character frequencies':
                export_data = _char_freq_option()
            elif input.download_type() == '% recovery':
                export_data = _percent_recovery_option()
            else:
                export_data = (None, None)

            # Create a temporary file for the download
            with tempfile.NamedTemporaryFile(prefix=export_data[1], suffix=f'.{download_format}', delete=False) as tmpfile:
                tmpfile.write(export_data[0].encode('utf-8'))
                tmpfile.flush()  # Ensure data is written to disk

                return tmpfile.name

        # catch other errors and display them
        except (FileNotFoundError, ValueError) as error:
            ui.notification_show(ui.tags.div(
                str(error),
                style="color: red; font-weight: bold;"
            ), duration=10)

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

    #### Handle everything in the analysis tab ####
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
    @reactive.event(input.analysis_plot_type_left)
    def update_additional_options_left():
        """
        Update UI for the left plot in the analysis tab
        """
        # ensure that it is switched back
        if input.analysis_plot_type_left() == 'Off':
            ui.remove_ui(selector="div:has(> #additional_analysis_options_left)")
            ui.remove_ui(selector="div:has(> #additional_analysis_options_left-label)")
            ui.remove_ui(selector="div:has(> #analysis_info_left)")
        if input.analysis_plot_type_left() == 'Pairwise identity':
            ui.insert_ui(
                ui.input_selectize(
                    'additional_analysis_options_left',
                    label='Options left',
                    choices={
                        'ghd': 'global hamming distance',
                        'lhd': 'local hamming distance',
                        'ged': 'gap excluded distance',
                        'gcd': 'gap compressed distance'
                    },
                    selected='ghd'
                ),
                selector='#analysis_plot_type_right-label',
                where='beforeBegin'
            )
            ui.insert_ui(
                ui.div(
                    ui.output_text_verbatim(
                        'analysis_info_left', placeholder=False
                    )
                ),
                selector='#analysis_plot_type_right-label',
                where='beforeBegin'
            )

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

    @render_widget
    def analysis_custom_heatmap():
        """
        Create the heatmap
        """
        aln = reactive.alignment.get()
        inputs = prepare_minimal_inputs(left_plot=True, window_size=True)
        try:
            if inputs['analysis_plot_type_left'] == 'Pairwise identity':
                return create_analysis_custom_heatmap(aln, inputs)
            else:
                return None
        except TypeError:
            return None

    @render_widget
    def analysis_char_freq_heatmap():
        """
        Create character frequency heatmap
        """
        aln = reactive.alignment.get()
        inputs = prepare_minimal_inputs(right_plot=True, window_size=True)
        try:
            if inputs['analysis_plot_type_right'] == 'Recovery':
                return create_recovery_heatmap(aln, inputs)
            elif inputs['analysis_plot_type_right'] == 'Character frequencies':
                return create_freq_heatmap(aln,  inputs)
            return None
        except TypeError:
            return None

