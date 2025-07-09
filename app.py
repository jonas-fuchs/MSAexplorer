"""
This contains the code to create the MSAexplorer shiny application
"""

# build-in
from pathlib import Path
# libs
from shiny import App
# app resource
from app_src.shiny_user_interface import shiny_ui
from app_src.shiny_server import server

# file paths for css and js
css_file = Path(__file__).resolve().parent/'app_src'/'www'/'css'/'styles.css'
js_file = Path(__file__).resolve().parent/'app_src'/'www'/'js'/'helper_functions.js'

# create the app
app = App(
    shiny_ui(
        css_file=css_file,
        js_file=js_file
    ),
    server,
    static_assets={'/img': Path(__file__).parent/'app_src'/'img'}
)
