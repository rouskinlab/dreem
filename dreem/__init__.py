__version__ = '0.0.5'

import warnings

from . import align, demultiplex, vector, cluster, aggregate

warnings.simplefilter(action='ignore', category=FutureWarning)
