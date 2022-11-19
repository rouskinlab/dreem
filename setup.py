from pkg_resources import Requirement
from setuptools import setup
import os, sys
from dreem import __version__
import sys

try:
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
except:
    with open('../requirements.txt') as f:
        requirements = f.read().splitlines()

PYTHON_VERSION = (3,10)

if sys.version_info < PYTHON_VERSION:
    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")


setup(
   name='dreem',
   version=__version__,
   license="MIT",
   description='Yves Martin\', Scott Grote\'s and Matty Allan\'s implementation of Prof. Silvi Rouskin\'s DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   long_description= 'TODO',
   packages=['dreem'],  #same as name
   package_dir={'dreem': 'dreem'},
   py_modules=[
         'dreem',
            'dreem.demutliplexing',
            'dreem.alignment',
            'dreem.vectoring',
            'dreem.clustering',
            'dreem.aggregate',
            'dreem.post_processing',
            'dreem.draw',
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
    entry_points = {
    'console_scripts' : [
        'dreem = dreem.run : run', 
        'dreem-demultiplexing = dreem.demultiplexing.run : run',  
        'dreem-alignment = dreem.alignment.run : run',
        'dreem-vectoring = dreem.vectoring.run : run',
        'dreem-clustering = dreem.clustering.run : run',
        'dreem-aggregate = dreem.aggregate.run : run',
        'dreem-post-processing = dreem.post_processing.run : run'
    ]
}
)