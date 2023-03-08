from setuptools import setup, find_packages
#from dreem import __version__
import sys

print("Python version: {}".format(sys.version))
requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())


#PYTHON_VERSION = (3,10)

#if sys.version_info < PYTHON_VERSION:
#    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
  long_description=readme,
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
   install_requires=requirements,
   entry_points = {'console_scripts' : ['dreem = dreem.main : cli']}
)
