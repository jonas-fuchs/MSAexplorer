"""
A script to provide minimal argument parsing for displaying the version or
launching the MSAexplorer app. Not needed for usage as a python module.

This module defines functions to parse command-line arguments and either display
the version information or start the MSAexplorer application. It primarily serves
as an entry point for launching the application or retrieving its version.
"""

import sys
import argparse
from msaexplorer import __version__


def parse_args(sysargs):
    """
    Minimal argument parser for displaying the version or launching the app.
    """
    parser = argparse.ArgumentParser(
        description='The MSAexplorer app is an interactive visualization tool designed for exploring multiple sequence alignments (MSAs).',
        usage='''\tmsaexplorer --run'''
    )

    parser.add_argument(
         '--run',
        action='store_true',
        help='Start the MSAexplorer app'
    )

    parser.add_argument(
        '--host',
        type=str,
        default='127.0.0.1',
        metavar='ip',
        help='The address that the app should listen on. Defaults to 127.0.0.1'
    )

    parser.add_argument(
        '--port',
        type=int,
        default=8080,
        metavar='port',
        help='The port that the app should listen on. Set to 0 to use a random port. Defaults to 8080.'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'MSAexplorer {__version__}'
    )

    if not sysargs:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args(sysargs)


def main(sysargs=sys.argv[1:]):
    args = parse_args(sysargs)

    if args.run:
        try:
            from shiny import run_app
            from shiny import App
            from app_src.shiny_user_interface import shiny_ui
            from app_src.shiny_server import server
            from importlib.resources import files
        except ImportError:
            sys.exit(
                "Please install the MSAexplorer front end app via 'pip install msaexplorer[app]' or 'pip install msaexplorer[app-plus]'."
            )

        css_path = files("app_src").joinpath("www/css/styles.css")
        js_path = files("app_src").joinpath("www/js/helper_functions.js")
        img_path = files("app_src").joinpath("www/img")
        # same code as in root/app.py
        run_app(
            app=App(
                shiny_ui(
                    css_file=css_path,
                    js_file=js_path
                ),
                server,
                static_assets={'/img': str(img_path)}
            ),
            port=args.port,
            host=args.host
        )


if __name__ == "__main__":
    main()
