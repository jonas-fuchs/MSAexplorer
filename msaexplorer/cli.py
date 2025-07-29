"""
A script to provide minimal argument parsing for displaying the version or
launching the MSAexplorer app. Not needed for usage as a python module.

This module defines functions to parse command-line arguments and either display
the version information or start the MSAexplorer application. It primarily serves
as an entry point for launching the application or retrieving its version.
"""

import sys
import argparse
from pathlib import Path
from shiny import run_app
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
        app_path = Path(__file__).parent.parent / "app_src" / "app.py"
        run_app(str(app_path))


if __name__ == "__main__":
    main()
