from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from dreem.util.cli import CountOption, CountOptionValue
from dreem.util.util import AMBIG_INT
from dreem.vector.mprofile import VectorReader, VectorTextTranslator



def run(out_dir: str,
        report: tuple[str],
        count: tuple[str] = (),
        frac: tuple[str] = ()):
    for rep in report:
        reader = VectorReader.from_report_file(rep)
        with open(f"test_vectors.txt", "wb") as f:
            f.write(VectorTextTranslator.blocktrans(reader.get_all_vectors()))

        for query in count:
            quint = CountOptionValue[CountOption(query).name]
            colors = [{"A": "red", "C": "blue", "G": "orange", "T": "green"}[reader.ref_seq.decode()[pos-1]]
                      for pos in reader.positions]
            counts = reader.count_query(quint)
            plt.bar(counts.index, counts, color=colors)
            plt.show()
        for query in frac:
            quint = CountOptionValue[CountOption(query).name]
            colors = [{"A": "red", "C": "blue", "G": "orange", "T": "green"}[reader.ref_seq.decode()[pos-1]]
                      for pos in reader.positions]
            counts = reader.count_query(quint)
            xcounts = reader.count_query(AMBIG_INT ^ quint)
            f = counts / (counts + xcounts)
            plt.bar(f.index, f, color=colors)
            plt.show()
