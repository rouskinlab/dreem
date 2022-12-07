"""Init file """
__version__= '0.0.1'

import sys, os

#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dreem'))

import dreem.demultiplex.main as demultiplexing
import dreem.align.main as alignment
import dreem.vector.main as vectoring
import dreem.cluster.main as clustering
import dreem.aggregate.main as aggregation
import dreem.draw.main as drawer


