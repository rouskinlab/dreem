from setuptools import setup, find_packages

from Cython.Build import cythonize
import Cython.Compiler.Options


# Generate an HTML report for each Cythonized module.
Cython.Compiler.Options.annotate = True


requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())

# PYTHON_VERSION = (3,10)

# if sys.version_info < PYTHON_VERSION:
#    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
    name="dreem",
    version='0.1.1',
    description="DREEM solves RNA structure ensembles using chemical probing data",
    long_description=readme,
    author="Silvi Rouskin Lab",
    author_email="silvi@hms.harvard.edu",
    url="https://github.com/rouskinlab/dreem",
    packages=find_packages(),
    package_dir={'dreem': 'dreem'},
    #py_modules=[
    #    'dreem',
    #    'dreem/demultiplex',
    #    'dreem/align',
    #    'dreem/vector',
    #    'dreem/cluster',
    #    'dreem/aggregate',
    #    'dreem/draw',
    #    'dreem/test',
    #    'dreem/util',
    #],
    ext_modules=cythonize("dreem/vector/vector.py", annotate=True),
    include_package_data=True,
    package_data={
        "dreem": ["test-data/vector-test-data/vectorize-read-test-data.csv"]
    },
    install_requires=requirements,
    entry_points={'console_scripts': ['dreem = dreem.main : cli']},
    python_requires=">=3.10, <4",
)
