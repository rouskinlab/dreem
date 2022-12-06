
import json
import os,sys
sys.path.append(os.getcwd())
import dreem.util as util
import dreem
import subprocess
import pytest

def run_notebook(notebook):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    try: 
        os.system('jupyter nbconvert --to notebook --execute '+ os.path.join(current_dir,'test_files_generators', notebook))
    except:
        raise 'Error running notebook {}'.format(notebook)
    try:
        os.remove(os.path.join(current_dir,'test_files_generators', notebook.split('/')[-1].split('.')[0]+'.nbconvert.ipynb'))
    except:
        raise 'Error while removing nbconvert file. Issue is probably that the notebook hasn\'t been run properly'

def test_make_demultiplexing_files():
    run_notebook('make_demultiplexing_files.ipynb')

def test_make_alignment_files():
    run_notebook('make_alignment_files.ipynb')

def test_make_vectoring_files():
    run_notebook('make_vectoring_files.ipynb')

def test_make_clustering_files():
    run_notebook('make_clustering_files.ipynb')

def test_make_aggregate_files():
    run_notebook('make_aggregate_files.ipynb')

# remove this test if you want to run the tests!!
def test_make_predicted_output_an_output():
    root = os.path.abspath(os.path.join(os.path.realpath(__file__),'..','test_files'))
    print(root)
    os.system('rm -r {}/output'.format(root))
    os.system('cp -fr {}/predicted_output {}/output'.format(root,root))
    assert os.path.exists(os.path.join(root,'output')), 'The output folder doesn\'t exist'