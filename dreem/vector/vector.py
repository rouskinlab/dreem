"""
Python interface to the compiled C extension module for vectoring.
"""

import os
import sys


try:
    # Check if the C vectoring module is already installed.
    from vectorc import *
except ModuleNotFoundError:
    # If not, then see if the project has been built.
    vector_dir = os.path.dirname(__file__)
    source_dir = os.path.dirname(vector_dir)
    project_dir = os.path.dirname(source_dir)
    build_dir = os.path.join(project_dir, "build")
    vector_so_modules = list()
    try:
        print(f"Searching {build_dir} for modules ...")
        # If the project has been built, then look for the module.
        for dirname, subdirs, files in os.walk(build_dir):
            print(dirname)
            vector_so_modules.extend(os.path.join(dirname, f) for f in files
                                     if (f.startswith("vectorc")
                                         and f.endswith(".so")))
    except FileNotFoundError:
        pass
    print("Found:", vector_so_modules)
    # Check if the module was found and if it is unique.
    if len(vector_so_modules) == 0:
        raise ModuleNotFoundError("No module named 'vectorc'")
    if len(vector_so_modules) > 1:
        raise ImportError("Multiple modules named 'vectorc'")
    # If so, then add it to the path.
    vector_so_module = vector_so_modules[0]
    sys.path.append(os.path.dirname(vector_so_module))
    # Attempt to import the module.
    from vectorc import *

import vectorc
print(f"Imported vectorc module in file: {vectorc.__file__}")

