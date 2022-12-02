
import os
import pandas as pd

#output_folder = os.path.join(OUTPUT_FOLDER, 'vectoring')

def test_done():
    assert False, 'Test isn\'t written yet'

def test_repo_exists():
    assert os.path.exists(os.path.join(output_folder, SAMPLE_NAME)), 'Directory {} does not exist'.format(os.path.join(output_folder, SAMPLE_NAME))

def index_starts_with_0(directory):
	for construct in CONSTRUCTS:
		bv_file = os.path.join(directory, SAMPLE_NAME, construct)+'.orc'
		assert os.path.exists(bv_file), 'Bitvector file does not exist for construct {}'.format(construct)
		df = pd.read_orc(bv_file)
		assert list(df.columns)[0][1] == '0', 'Index does not start with 0'

def test_index_starts_with_0():
	#directory= os.path.join(output_folder, SAMPLE_NAME)
	#index_starts_with_0(directory)
	pass
