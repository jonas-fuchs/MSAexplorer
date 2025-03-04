#TODO: add function to run shiny
import subprocess
from pathlib import Path

app = Path(__file__).parent / 'app.py'

def main():
    subprocess.run(f'shiny run {app}', capture_output=True, check=True)