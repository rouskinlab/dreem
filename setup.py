from setuptools import setup, find_packages
from src import __version__
import sys

requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())


PYTHON_VERSION = (3,10)

if sys.version_info < PYTHON_VERSION:
    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
   name='src',
   version=__version__,
   license="MIT",
   description=("Implementation of Prof Silvi Rouskin's DREEM algorithm "
                "by Yves Martin, Scott Grote, and Matty Allan"),
   author="Yves Martin des Taillades, Scott Grote, and Matty Allan",
   author_email='yves@martin.yt',
   long_description=readme,
   url='https://github.com/rouskinlab/dreem',
   packages=find_packages(),
   package_dir={'src': 'src'},
   py_modules=[
        'src',
        'src/demultiplex',
        'src/align',
        'src/vector',
        'src/cluster',
        'src/aggregate',
        'src/draw',
        'src/test',
        'src/util',
   ],
   include_package_data=True,
   install_requires=requirements,
   entry_points = {'console_scripts' : ['dreem = src.main : cli']}
)
