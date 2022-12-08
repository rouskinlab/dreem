
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import subprocess, os

def test_remove_files():
    subprocess.run(['rm', '-fr', test_files_dir])