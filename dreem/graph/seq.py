from abc import ABC, abstractmethod
from functools import cache
from itertools import chain
import os
from pathlib import Path

from click import command
import pandas as pd
from plotly import graph_objects as go

from .base import (PRECISION, find_tables, GraphWriter, CartesianGraph,
                   OneTableSeqGraph, OneSampGraph, OneTableSectGraph)
from .color import RelColorMap, SeqColorMap
from ..core import docdef, path
from ..core.cli import (opt_table, opt_fields, opt_stack, opt_yfrac,
                        opt_csv, opt_html, opt_pdf, opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.seq import BASES
from ..table.base import CountTable
from ..table.load import (POS_FIELD, TableLoader, RelPosTableLoader,
                          MaskPosTableLoader, ClusterPosTableLoader)

# Number of digits to which to round decimals.

params = [
    opt_table,
    opt_fields,
    opt_yfrac,
    opt_stack,
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
        fields: str,
        yfrac: bool,
        stack: bool, *,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    writers = list(map(SeqGraphWriter, find_tables(table)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(fields=fields, count=yfrac,
                                            stack=stack, csv=csv,
                                            html=html, pdf=pdf))))


class SeqGraphWriter(GraphWriter):

    def iter(self, fields: str, count: bool, stack: bool):
        if isinstance(self.table, RelPosTableLoader):
            if stack:
                if count:
                    yield RelCountStackedPosBarGraph(table=self.table,
                                                     codes=fields)
                else:
                    yield RelFracStackedPosBarGraph(table=self.table,
                                                    codes=fields)
            else:
                for field in fields:
                    if count:
                        yield RelCountSerialPosBarGraph(table=self.table,
                                                        codes=field)
                    else:
                        yield RelFracSerialPosBarGraph(table=self.table,
                                                       codes=field)
        elif isinstance(self.table, MaskPosTableLoader):
            if stack:
                if count:
                    yield MaskCountStackedPosBarGraph(table=self.table,
                                                      codes=fields)
                else:
                    yield MaskFracStackedPosBarGraph(table=self.table,
                                                     codes=fields)
            else:
                for field in fields:
                    if count:
                        yield MaskCountSerialPosBarGraph(table=self.table,
                                                         codes=field)
                    else:
                        yield MaskFracSerialPosBarGraph(table=self.table,
                                                        codes=field)
        elif isinstance(self.table, ClusterPosTableLoader):
            for cluster in self.table.cluster_names:
                yield ClustPosBarGraph(table=self.table, cluster=cluster)


class SeqGraph(CartesianGraph, OneTableSeqGraph, OneSampGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *args,
                 table: (RelPosTableLoader
                         | MaskPosTableLoader
                         | ClusterPosTableLoader),
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table

    @classmethod
    @abstractmethod
    def get_source(cls):
        """ Step from which the data came. """
        return ""

    @classmethod
    def get_table_type(cls):
        return RelPosTableLoader, MaskPosTableLoader, ClusterPosTableLoader

    @classmethod
    def get_data_type(cls):
        return pd.Series

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @classmethod
    def get_xattr(cls):
        return POS_FIELD


class SeriesSeqGraph(SeqGraph, ABC):
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


class StackedSeqGraph(SeqGraph, ABC):
    """ Bar graph of a table with multiple types of bit fields, with one
    bar for each position in the sequence. """

    '''
    def __init__(self, *args, codes:  **kwargs):
        super().__init__(*args, **kwargs)
        self.

    @classmethod
    @abstractmethod
    def _validate_codes(cls, codes: str):
        """ Confirm that the codes are valid. """

    @codes.setter
    def codes(self, codes: str):
        self._validate_codes(codes)
        self._codes = codes
    
    '''

    @classmethod
    @abstractmethod
    def get_yattr(cls):
        return ""

    @classmethod
    @abstractmethod
    def get_table_field(cls, _: CountTable | TableLoader, __: str):
        """ Load one field from the table. """
        return pd.Series()

    @property
    def title(self):
        fields = '/'.join(sorted(CountTable.FIELD_CODES[c] for c in self.codes))
        return (f"{self.sample}: {self.get_source()} {self.get_yattr()}s of "
                f"{fields} bases per position in {self.ref}")

    @property
    def graph_filename(self):
        sort_codes = "".join(sorted(self.codes))
        by = "stacked" if self.is_stacked() else "serial"
        fname = f"{self.get_source()}_{self.get_yattr()}_{sort_codes}_{by}"
        return fname.lower()


class CountFieldPosBarGraph(StackedSeqGraph, ABC):
    """ FieldPosBarGraph where each bar represents a count of reads. """

    @classmethod
    def get_yattr(cls):
        return "Count"

    @classmethod
    def get_table_field(cls, table: CountTable | TableLoader, code: str):
        return table.get_field_count(code)


class FracFieldPosBarGraph(StackedSeqGraph, ABC):
    """ FieldPosBarGraph where each bar represents a fraction of reads. """

    @classmethod
    def get_yattr(cls):
        return "Fraction"

    @classmethod
    def get_table_field(cls, table: CountTable | TableLoader, code: str):
        return table.get_field_frac(code).round(PRECISION)


class SeriesFieldPosBarGraph(StackedSeqGraph, SeriesSeqGraph, ABC):
    """ Bar graph wherein each bar represents a base in a sequence. """

    def __init__(self, *args, codes: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.codes = codes

    @classmethod
    def is_stacked(cls):
        return False

    @classmethod
    def _validate_codes(cls, codes: str):
        if len(codes) != 1:
            raise ValueError(f"Expected 1 field code, but got {len(codes)}")

    @property
    def code(self):
        """ The code of the field to graph. """
        return self.codes[0]

    def _get_data(self) -> pd.Series:
        return self.get_table_field(self.table, self.code)


class StackedFieldPosBarGraph(StackedSeqGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    outcomes for a base in a sequence. """

    def __init__(self, *args, codes: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.codes = codes

    @classmethod
    def is_stacked(cls):
        return True

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    @classmethod
    def get_data_type(cls):
        return pd.DataFrame

    @classmethod
    def _validate_codes(cls, codes: str):
        if len(codes) == 0:
            raise ValueError("Expected 1 or more field codes, but got 0")

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


class RelPosBarGraph(StackedSeqGraph, ABC):
    """ Bar graph of related data from one sample at each position in a
    sequence. """

    @classmethod
    def get_source(cls):
        return "Related"


class MaskPosBarGraph(StackedSeqGraph, SectSeqGraph, ABC):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """

    @classmethod
    def get_source(cls):
        return "Masked"


class RelFracSerialPosBarGraph(SeriesFieldPosBarGraph,
                               FracFieldPosBarGraph,
                               RelPosBarGraph):
    """ """


class RelFracStackedPosBarGraph(StackedFieldPosBarGraph,
                                FracFieldPosBarGraph,
                                RelPosBarGraph):
    """ """


class RelCountSerialPosBarGraph(SeriesFieldPosBarGraph,
                                CountFieldPosBarGraph,
                                RelPosBarGraph):
    """ """


class RelCountStackedPosBarGraph(StackedFieldPosBarGraph,
                                 CountFieldPosBarGraph,
                                 RelPosBarGraph):
    """ """


class MaskFracSerialPosBarGraph(SeriesFieldPosBarGraph,
                                FracFieldPosBarGraph,
                                MaskPosBarGraph):
    """ """


class MaskFracStackedPosBarGraph(StackedFieldPosBarGraph,
                                 FracFieldPosBarGraph,
                                 MaskPosBarGraph):
    """ """


class MaskCountSerialPosBarGraph(SeriesFieldPosBarGraph,
                                 CountFieldPosBarGraph,
                                 MaskPosBarGraph):
    """ """


class MaskCountStackedPosBarGraph(StackedFieldPosBarGraph,
                                  CountFieldPosBarGraph,
                                  MaskPosBarGraph):
    """ """


class ClustPosBarGraph(SeriesSeqGraph, SectSeqGraph):
    """ Bar graph of a table of per-cluster mutation rates. """

    def __init__(self, *args, cluster: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.cluster = cluster

    @classmethod
    def get_source(cls):
        return "Clustered"

    @classmethod
    def get_yattr(cls):
        return "Mutation Rate"

    @property
    def title(self):
        return (f"{self.sample}: {self.get_source()} {self.get_yattr()}s of "
                f"bases per position in {self.ref}:{self.sect}, {self.cluster}")

    def _get_data(self):
        return self.table.data[self.cluster]

    @property
    def graph_filename(self):
        return path.fill_whitespace(self.cluster).lower()
