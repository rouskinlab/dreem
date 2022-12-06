
import json
import os,sys
# add the current directory to the path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files_generators'))
import dreem.util as util

from test_files_generators import aggregate, clustering, demultiplexing, vectoring, alignment, main

def test_make_demultiplexing_files():
    demultiplexing.make_files()

def test_make_alignment_files():
    alignment.make_files()

def test_make_vectoring_files():
    vectoring.make_files()

def test_make_clustering_files():
    clustering.make_files()

def test_make_aggregate_files():
    aggregate.make_files()
    
def test_make_main_files():
    main.make_files()

# remove this test if you want to run the tests!!
def test_make_predicted_output_an_output():
    root = os.path.abspath(os.path.join(os.path.realpath(__file__),'..','test_files'))
    os.system('rm -r {}/output'.format(root))
    os.mkdir(os.path.join(root,'output'))
    for module in ['alignment','vectoring','clustering','aggregate', 'main']:
        os.system('rm -r {}/output/{}'.format(root, module))
        os.system('cp -fr {}/predicted_output/{} {}/output/{}'.format(root, module, root, module))
        assert os.path.exists(os.path.join(root,'predicted_output',module)), 'The predicted output folder doesn\'t exist'