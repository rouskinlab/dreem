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
   name='dreem_nap',
   version=__version__,
   license="MIT",
   description='Yves Martin and Scott Grote\' implementation of Prof. Silvi Rouskin\'s DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   long_description= 'TODO',
   packages=['dreem'],  #same as name
   package_dir={'dreem': 'dreem'},
   py_modules=[
       'dreem/draw',
       'dreem/draw/util', 
       'dreem/draw/manipulator',
       'dreem/draw/study',
       'dreem/draw/plotter',
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
   
)