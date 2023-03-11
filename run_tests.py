import dreem
import os

test_pipeline_path = os.path.join(os.path.dirname(dreem.__file__), 'test', 'test_pipeline.py')

os.system("pytest " + test_pipeline_path)
