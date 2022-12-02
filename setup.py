from setuptools import setup, find_packages
from dreem import __version__
import sys

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


PYTHON_VERSION = (3,10)

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
            'dreem/demultiplexing',
            'dreem/alignment',
            'dreem/vectoring',
            'dreem/clustering',
            'dreem/aggregate',
            'dreem/draw',
            'dreem/util',
            'test',
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
    ]
}
)