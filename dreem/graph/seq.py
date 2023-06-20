from abc import ABC, abstractmethod
from functools import cache
from itertools import chain
from logging import getLogger
import os
from pathlib import Path

from click import command
import pandas as pd
from plotly import graph_objects as go

from .base import (PRECISION, find_tables, GraphWriter, CartesianGraph,
                   OneTableSeqGraph, OneSampGraph, OneTableSectGraph)
from .color import RelColorMap, SeqColorMap
from ..core import docdef
from ..core.cli import (opt_table, opt_rels, opt_stacks, opt_yfrac,
                        opt_csv, opt_html, opt_pdf, opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.seq import BASES
from ..table.base import Table
from ..table.load import (POS_FIELD, TableLoader, RelPosTableLoader,
                          MaskPosTableLoader, ClusterPosTableLoader)

logger = getLogger(__name__)

# Number of digits to which to round decimals.

params = [
    opt_table,
    opt_rels,
    opt_stacks,
    opt_yfrac,
    opt_csv,
    opt_html,
    opt_pdf,
    opt_max_procs,
    opt_parallel,
]


@command(__name__.split(os.path.extsep)[-1], params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positional attributes. """
    return run(*args, **kwargs)


@docdef.auto()
def run(table: tuple[str, ...],
        rels: str,
        stacks: bool,
        yfrac: bool, *,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    writers = list(map(SeqGraphWriter, find_tables(table)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels=rels, stacks=stacks,
                                            yfrac=yfrac, csv=csv,
                                            html=html, pdf=pdf))))


class SeqGraphWriter(GraphWriter):

    def iter(self, rels: str, stacks: str, yfrac: bool):
        if isinstance(self.table, RelPosTableLoader):
            for rel in rels:
                yield SerialRelSeqGraph(table=self.table,
                                        codes=rel,
                                        yfrac=yfrac)
            if stacks:
                yield StackedRelSeqGraph(table=self.table,
                                         codes=stacks,
                                         yfrac=yfrac)
        elif isinstance(self.table, MaskPosTableLoader):
            for rel in rels:
                yield SerialMaskSeqGraph(table=self.table,
                                         codes=rel,
                                         yfrac=yfrac)
            if stacks:
                yield StackedMaskSeqGraph(table=self.table,
                                          codes=stacks,
                                          yfrac=yfrac)
        else:
            logger.error(f"{self.__class__.__name__} cannot graph {self.table}")


class SeqGraph(CartesianGraph, OneTableSeqGraph, OneSampGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *args,
                 table: (RelPosTableLoader
                         | MaskPosTableLoader
                         | ClusterPosTableLoader),
                 codes: str,
                 yfrac: bool,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table
        self.codes = codes
        self.yfrac = yfrac

    @property
    def source(self):
        if isinstance(self.table, RelPosTableLoader):
            return "Related"
        if isinstance(self.table, MaskPosTableLoader):
            return "Masked"
        if isinstance(self.table, ClusterPosTableLoader):
            return "Clustered"
        raise TypeError(
            f"Invalid table type for {self}: {type(self.table).__name__}")

    @classmethod
    def get_data_type(cls):
        return pd.Series

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    def get_xattr(self):
        return POS_FIELD

    def get_yattr(self):
        return "Fraction" if self.yfrac else "Count"

    @property
    def title(self):
        fields = '/'.join(sorted(Table.FIELD_CODES[c] for c in self.codes))
        return (f"{self.get_yattr()}s of {fields} bases in {self.source} reads "
                f"from {self.sample} per position in {self.ref}")

    @property
    @abstractmethod
    def sort_codes(self):
        return "".join(sorted(self.codes))

    @property
    def graph_filename(self):
        return f"{self.source}_{self.sort_codes}_{self.get_yattr()}".lower()

    def get_table_field(self, table: Table | TableLoader, code: str):
        return (table.get_field_frac(code).round(PRECISION) if self.yfrac
                else table.get_field_count(code))


class SerialSeqGraph(SeqGraph, ABC):
    """ Bar graph with a single series of data. """

    def get_traces(self):
        traces = list()
        # Construct a trace for each type of base.
        for base in BASES:
            # Find the non-missing value at every base of that type.
            vals = self.data.loc[self.seqarr == base].dropna()
            # Check if there are any values to graph.
            if vals.size > 0:
                # Define the text shown on hovering over a bar.
                hovertext = [f"{chr(base)}{x}: {y}" for x, y in vals.items()]
                # Create a trace comprising all bars for this base type.
                traces.append(go.Bar(name=chr(base), x=vals.index, y=vals,
                                     marker_color=self.cmap[base],
                                     hovertext=hovertext,
                                     hoverinfo="text"))
        return traces

    @property
    def code(self):
        """ The code of the field to graph. """
        if len(self.codes) != 1:
            raise ValueError(f"Expected 1 code, but got {len(self.codes)}")
        return self.codes[0]

    @property
    def sort_codes(self):
        return "-".join(["serial", super().sort_codes])

    def _get_data(self) -> pd.Series:
        return self.get_table_field(self.table, self.code)


class StackedSeqGraph(SeqGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    outcomes for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    @classmethod
    def get_data_type(cls):
        return pd.DataFrame

    @property
    def sort_codes(self):
        return "-".join(["stacked", super().sort_codes])

    def _get_data(self):
        data = dict()
        for code in self.codes:
            series = self.get_table_field(self.table, code)
            data[series.name] = series
        return pd.DataFrame.from_dict(data)

    def get_traces(self):
        traces = list()
        # Construct a trace for each field.
        for field, vals in self.data.items():
            # Define the text shown on hovering over a bar.
            hovertext = [f"{chr(base)}{x} {field}: {y}"
                         for base, (x, y) in zip(self.seq, vals.items(),
                                                 strict=True)]
            # Create a trace comprising all bars for this field.
            traces.append(go.Bar(name=field, x=vals.index, y=vals,
                                 marker_color=self.cmap[field],
                                 hovertext=hovertext,
                                 hoverinfo="text"))
        return traces

    @cache
    def get_figure(self):
        fig = super().get_figure()
        # Stack the bars at each position.
        fig.update_layout(barmode="stack")
        return fig


class SectSeqGraph(SeqGraph, OneTableSectGraph, ABC):
    """ Bar graph of the positions in a section. """

    @property
    def title(self):
        return f"{super().title}:{self.sect}"


class RelSeqGraph(SeqGraph, ABC):
    """ """


class MaskSeqGraph(SectSeqGraph, ABC):
    """ """


class ClustSeqGraph(SectSeqGraph, ABC):
    """ """


class SerialRelSeqGraph(SerialSeqGraph, RelSeqGraph):
    """ Bar graph of related data from one sample at each position in a
    sequence. """


class StackedRelSeqGraph(StackedSeqGraph, RelSeqGraph):
    """ Bar graph of related data from one sample at each position in a
    sequence. """


class SerialMaskSeqGraph(SerialSeqGraph, MaskSeqGraph, ABC):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """


class StackedMaskSeqGraph(StackedSeqGraph, MaskSeqGraph, ABC):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """
