"""Init file """
__version__= '0.0.5'


import dreem.demultiplex.main as demultiplexing
import dreem.align.main as alignment
import dreem.vector.main as vectoring
import dreem.cluster.main as clustering
import dreem.aggregate.main as aggregation
#import dreem.draw.main as drawer

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
