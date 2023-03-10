__version__ = '0.1.0'

from multiprocessing import current_process
import random
import warnings

about = f"""
########################################################################

DREEM: Detection of RNA folding Ensembles using Expectation–Maximization
                         Version: {__version__}                         

          DREEM was developed in the lab of Dr. Silvi Rouskin           
      Assistant Professor of Microbiology, Harvard Medical School

The original version of DREEM is described in the following publication:
Tomezsko, P.J., Corbin, V.D.A., Gupta, P. et al. Determination of RNA
structural diversity and its role in HIV-1 RNA splicing. Nature 582,
438–442 (2020). https://doi.org/10.1038/s41586-020-2253-5

This version was written by the following affiliates of the Rouskin lab:
- Yves Martin des Taillades (Stanford University)
- Scott Grote (Harvard Medical School)
- Matthew "Matty" Allan (MIT, Harvard Medical School)
- Albéric de Lajarte (Stanford University)

########################################################################
"""

loading_messages = [
    "Please wait while DREEM loads",
    "DREEM is loading ... please wait",
    "Loading DREEM ... please wait"
]

if current_process().name == "MainProcess":
    # Print the introductory message only if this is the main process,
    # not a child process created by the multiprocessing module.
    message = loading_messages[random.randint(0, len(loading_messages) - 1)]
    print(f"{message.strip()}\n{about.strip()}\n")

from . import align, demultiplex, vector, cluster, aggregate, draw

from .main import run

warnings.simplefilter(action='ignore', category=FutureWarning)
