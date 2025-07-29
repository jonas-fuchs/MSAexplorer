"""
This contains the code to create the MSAexplorer shiny application with shinylive. Not needed for usage as a python module.
"""

# build-in
from importlib.resources import files

# libs
from shiny import App
# app resource
from app_src.shiny_user_interface import shiny_ui
from app_src.shiny_server import server

css_path = files("app_src").joinpath("www/css/styles.css")
js_path = files("app_src").joinpath("www/js/helper_functions.js")
img_path = files("app_src").joinpath("www/img")

# create the app
app = App(
    shiny_ui(
        css_file=css_path,
        js_file=js_path
    ),
    server,
    static_assets={'/img': str(img_path)}
)
