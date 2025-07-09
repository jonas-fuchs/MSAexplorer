"""
This module contains the resources to generate the main ui.
"""

# import libs
from shiny import ui
from shinywidgets import output_widget
from matplotlib import colormaps

def shiny_ui(css_file, js_file):
    """
    main UI generation
    """
    app_ui = ui.page_fluid(
        # include css and js
        ui.include_css(css_file),
        ui.include_js(js_file),
        # Custom sidebar
        ui.div(id="overlay-bg", onclick="toggleSidebar()"),
        _custom_sidebar(),
        # Main side
        ui.navset_bar(
            _upload_tab(),
            _plot_tab(),
            _analysis_tab(),
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

    return app_ui


def _custom_sidebar():
    """
    custom sidebar that overlays all
    """

    return ui.div(
    ui.h6(
        ui.img(src='img/settings.svg', height='20px'),
        ' SETTINGS', style="margin-top: 15px;"),
            ui.hr(),
            ui.h6('Statistic settings (entropy, similarity etc)'),
            ui.layout_columns(
                ui.input_numeric('rolling_avg', 'Rolling average', value=1, min=1),
                ui.div(
                    ui.HTML(
                        """
                        <div style="display: flex; flex-direction: column; align-items: flex-start;">
                            <label for="stat_color" style="font-size: 14px; margin-bottom: 5px;">Color for stat plot:</label>
                            <input type="color" id="stat_color" value="#808080" onchange="updateColor(this.value)" 
                                   style="width: 35px; height: 35px; padding: 0; border: 1px; margin-bottom: 15px;">
                        </div>
                        """
                    )
                ),
                ui.input_selectize('logo_coloring', 'Coloring sequence logo', ['standard'], selected='standard')
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
                    ui.input_switch('seq_names', 'show names', value=False),

                ),
                ui.column(
                    4,
                ui.input_selectize('reference', 'Reference', ['first', 'consensus'], selected='first'),
                    ui.div(
                        ui.HTML(
                            """
                            <div style="display: flex; flex-direction: column; align-items: flex-start;">
                                <label for="reference_color" style="font-size: 14px; margin-bottom: 5px;">Reference color:</label>
                                <input type="color" id="reference_color" value="#4682B4" onchange="updateColor(this.value)" 
                                       style="width: 35px; height: 35px; padding: 0; border: 1px; margin-bottom: 15px;">
                            </div>
                            """
                        )
                    ),
                    ui.input_switch('show_sequence', 'show sequence for differences', value=True),
                    ui.input_switch('show_sequence_all', 'show all sequences', value=False),
                ),
                ui.column(
                4,

                    ui.input_selectize('matrix', 'Matrix', ['None']),
                    ui.input_selectize('matrix_color_mapping', 'Colormap Similarity', choices=list(colormaps.keys()), selected='PuBu_r'),
                    ui.input_selectize('identity_coloring', 'Identity coloring', ['None', 'standard'], selected='None'),
                )
            ),
            ui.hr(),
            ui.h6('Annotation settings'),
            ui.input_switch('show_legend_third_plot', 'Legend', value=True),
            ui.row(
                ui.column(
                4,
                ui.h6('SNP plot'),
                    ui.input_numeric('head_size', 'Head size', value=3, min=1),
                    ui.input_numeric('stem_size', 'Stem length', value=1, min=1),
                    ui.input_selectize('snp_coloring', 'Coloring', ['standard'], selected='standard'),
                    id='snp_column'
                ),
                ui.column(
                4,
                ui.h6('Annotation plot'),
                    ui.input_selectize('feature_display', 'Feature', ['None']),
                    ui.div(
                        ui.HTML(
                            """
                            <div style="display: flex; flex-direction: column; align-items: flex-start;">
                                <label for="feature_color" style="font-size: 14px; margin-bottom: 5px;">Feature color:</label>
                                <input type="color" id="feature_color" value="#808080" onchange="updateColor(this.value)" 
                                       style="width: 35px; height: 35px; padding: 0; border: 1px; margin-bottom: 15px;">
                            </div>
                            """
                        )
                    ),
                    ui.input_numeric('strand_marker_size', 'Strand marker size', value=5, min=1, max=20)
                ),
            ),
            id='overlay-sidebar',
        )


def _upload_tab():
    """
    creates the upload panel
    """
    return ui.nav_panel(
        ' UPLOAD/DOWNLOAD',
        ui.div(
            ui.layout_columns(
                ui.card(
                    ui.card_header(ui.h6('Upload files:')),
                    ui.layout_columns(
                        ui.input_file('alignment_file', 'Multiple sequence alignment (*.fasta):', multiple=False,
                                      accept=['.fa', '.fasta', '.aln']),
                        ui.input_file('annotation_file', 'Optional annotation file (*.gff, *.bed, *.gb):', multiple=False,
                                      accept=['.gff', '.gff3', '.bed', '.gb']),
                    )
                ),
                ui.card(
                    ui.card_header(ui.h6('Download files:'),
                                   ui.popover(
                                       ui.span(
                                           ui.HTML(
                                               '<img src="img/gear.svg" alt="settings" style="height:16px; width:16px; position:absolute; top: 10px; right: 7px;">')
                                       ),
                                       ui.input_selectize('download_format', label='Format:', choices=[]),
                                )
                            ),
                    ui.input_selectize('download_type', label='Choose:', choices=['SNPs', 'consensus']),
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
                "MSAexplorer is an interactive visualization tool designed for exploring multiple sequence alignments (MSAs)."
            ),
            ui.p(
                ui.a(
                    "🔗 Full API documentation", href="https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer.html",
                      target="_blank"
                )
            ),
            ui.p(
                ui.a(
                    "🔗 Contribute on GitHub", href="https://github.com/jonas-fuchs/MSAexplorer",
                    target="_blank"
                )
            ),
            class_="about-card"
        ),
        icon=ui.HTML('<img src="img/upload.svg" alt="Upload Icon" style="height: 1em; vertical-align: middle">')
    ),


def _plot_tab():
    """
    creates the plot panel
    """
    return ui.nav_panel(
        ' PLOT',
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_slider('increase_height', 'Plot height', min=0.5, max=10, step=0.5, value=1),
                ui.input_selectize('stat_type', ui.h6('First plot'), ['Off'], selected='Off'),
                ui.input_numeric('plot_1_size', 'Plot fraction', 1, min=1, max=200),
                ui.input_selectize('alignment_type', ui.h6('Second plot'),['Off', 'identity', 'similarity'], selected='identity'),
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
            ui.row(
                ui.column(2, ui.input_numeric("zoom_start", "Start", 0, min=0, max=1000)),
                ui.column(8, ui.input_slider("zoom_range", "Zoom", min=0, max=1000, value=(0, 1000), step=1,
                                             width="100%")),
                ui.column(2, ui.input_numeric("zoom_end", "End", 1000, min=0, max=1000))
            ),
            ui.output_plot('msa_plot', height='100vh', width='92vw'),
            fillable=False
        ),
        icon=ui.HTML('<img src="img/chart.svg" alt="Chart Icon" style="height: 1em; vertical-align: middle;">')
    )


def _analysis_tab():
    """
    creates the analysis panel
    """
    return ui.nav_panel(
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
                ui.input_selectize('analysis_plot_type_left', ui.h6('Left plot'), ['Off', 'Pairwise identity'], selected='Off'),
                ui.input_selectize('analysis_plot_type_right', ui.h6('Right plot'), ['Off', 'Character frequencies', 'Recovery'], selected='Off'),
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
    )

