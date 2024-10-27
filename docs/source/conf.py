# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
#from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath('../../dreem/dreem'))



sys.path.append(os.path.abspath("..")+'/..')
#sys.path.insert(os.path.abspath('../..')) 	
sys.path.insert(0, os.path.abspath("../../"))
 
# Fix matplotlib non import
MOCK_MODULES = ['yaml']
with open('../../dreem/requirements.txt') as f:
    for line in f:
        MOCK_MODULES.append(line.strip().split('=')[0])
        
#for mod_name in MOCK_MODULES:
#    sys.modules[mod_name] = MagicMock()
    

# -- Project information -----------------------------------------------------

project = 'dreem'
copyright = '2023, Rouskin Lab'
author = 'Yves Martin des Taillades, Scott Grote, Matthew Allan, Alberic de Lajarte'

# The full version, including alpha/beta/rc tags
release = '28.02.2023'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

source_suffix = ['.rst', '.md']

# -- Options for HTML output
html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

autosectionlabel_prefix_document = True

# Options for build
# fail_on_warning = True
# autodoc_mock_imports = MOCK_MODULES
nitpick_ignore = [('py:class', 'type')]

# Generate the plots for the gallery
sys.path.append(os.path.abspath(""))
from plots import gallery_generator
gallery_generator.main()
