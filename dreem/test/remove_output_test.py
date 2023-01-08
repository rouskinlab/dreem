
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import subprocess, os

def test_remove_files():
    os.system(' '.join(['rm', '-fr', test_files_dir]))
    
def remove_module_files(module):
    for folder in ['input', 'output', 'expected_output']:
        os.system('rm -rf {}'.format(os.path.join(test_files_dir, folder, module)))
        os.makedirs(os.path.join(test_files_dir, folder, module), exist_ok=True)