from setuptools import setup, find_packages
from dreem import __version__
import sys

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


PYTHON_VERSION = (3,11)

if sys.version_info < PYTHON_VERSION:
    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
   name='dreem',
   version=__version__,
   license="MIT",
   description='Yves Martin\', Scott Grote\'s and Matty Allan\'s implementation of Prof. Silvi Rouskin\'s DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   long_description= readme,
   url='https://github.com/yvesmartindestaillades/dreem',
   packages=find_packages(),
   package_dir={'dreem': 'dreem'},
   py_modules=[
         'dreem',
            'dreem/demultiplex',
            'dreem/align',
            'dreem/vector',
            'dreem/cluster',
            'dreem/aggregate',
            'dreem/draw',
            'dreem/test',
            'dreem/util',
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
    entry_points = {
    'console_scripts' : [
        'dreem = dreem.cli : cli', 
        'dreem-demultiplexing = dreem.demultiplexing.cli : cli',  
        'dreem-alignment = dreem.alignment.cli : cli',
        'dreem-vectoring = dreem.vectoring.cli : cli',
        'dreem-clustering = dreem.clustering.cli : cli',
        'dreem-aggregate = dreem.aggregate.cli : cli',
    ]
}
)