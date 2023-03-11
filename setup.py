from setuptools import setup, find_packages
import sys

from Cython.Build import cythonize


requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())


PYTHON_VERSION = (3,10)

if sys.version_info < PYTHON_VERSION:
    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
    name='dreem',
    version="0.1.0",
    license="MIT",
    description=("Implementation of Prof Silvi Rouskin's DREEM algorithm "
                 "by Yves Martin, Scott Grote, and Matty Allan"),
    author="Yves Martin des Taillades, Scott Grote, and Matty Allan",
    author_email='yves@martin.yt',
    long_description=readme,
    url='https://github.com/rouskinlab/dreem',
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
    ext_modules=cythonize("dreem/vector/vector.pyx"),
    include_package_data=True,
    install_requires=requirements,
    entry_points = {'console_scripts' : ['dreem = dreem.main : cli']}
)
