
import json
import os,sys
sys.path.append(os.getcwd())
import dreem.util as util

    
def run_script(script):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    try: 
        os.system('python3 {}'.format(os.path.join(current_dir,'test_files_generators', script)))
    except:
        raise 'Error running script {}'.format(script)

def test_make_demultiplexing_files():
    run_script('make_demultiplexing_files.py')

def test_make_alignment_files():
    run_script('make_alignment_files.py')

def test_make_vectoring_files():
    run_script('make_vectoring_files.py')

def test_make_clustering_files():
    run_script('make_clustering_files.py')

def test_make_aggregate_files():
    run_script('make_aggregate_files.py')

# remove this test if you want to run the tests!!
def test_make_predicted_output_an_output():
    root = os.path.abspath(os.path.join(os.path.realpath(__file__),'..','test_files'))
    print(root)
    os.system('rm -r {}/output'.format(root))
    os.system('cp -fr {}/predicted_output {}/output'.format(root,root))
    assert os.path.exists(os.path.join(root,'output')), 'The output folder doesn\'t exist'