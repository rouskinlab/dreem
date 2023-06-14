from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
import os
from pathlib import Path
from typing import Iterable

from click import command
import pandas as pd
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from .base import find_tables, OneSampleRefGraph, OneSampleSectGraph
from .color import RelColorMap, SeqColorMap
from ..cluster.load import format_names_ks, parse_names
from ..core import docdef, path
from ..core.cli import (opt_table, opt_fields, opt_stack, opt_count,
                        opt_html, opt_pdf, opt_max_procs, opt_parallel)
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.seq import BASES
from ..table.base import CountTable, SEQ_FIELD
from ..table.load import (POS_FIELD, load, TableLoader, ClusterPosTableLoader,
                          MaskPosTableLoader, RelPosTableLoader)

# Number of digits to which to round decimals.
PRECISION = 6

params = [
    opt_table,
    opt_fields,
    opt_count,
    opt_stack,
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
        count: bool,
        stack: bool, *,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    return list(chain(*dispatch(graph_table, max_procs, parallel,
                                pass_n_procs=False,
                                args=as_list_of_tuples(find_tables(table)),
                                kwargs=dict(fields=fields, count=count,
                                            stack=stack, html=html, pdf=pdf))))


def iter_graphs(table: TableLoader, fields: str, count: bool, stack: bool):
    if isinstance(table, RelPosTableLoader):
        if stack:
            if count:
                yield StackCountRelatePosBarGraph.from_table(table, fields)
            else:
                yield StackFractionRelatePosBarGraph.from_table(table, fields)
        else:
            for field in fields:
                if count:
                    yield OneCountRelatePosBarGraph.from_table(table, field)
                else:
                    yield OneFractionRelatePosBarGraph.from_table(table, field)
    elif isinstance(table, MaskPosTableLoader):
        if stack:
            if count:
                yield StackCountMaskPosBarGraph.from_table(table, fields)
            else:
                yield StackFractionMaskPosBarGraph.from_table(table, fields)
        else:
            for field in fields:
                if count:
                    yield OneCountMaskPosBarGraph.from_table(table, field)
                else:
                    yield OneFractionMaskPosBarGraph.from_table(table, field)
    elif isinstance(table, ClusterPosTableLoader):
        if stack:
            yield ClustPosBarGraph.from_table(table)
        else:
            for k in table.ks:
                yield ClustPosBarGraph.from_table(table, [k])


def graph_table(file: Path,
                fields: str, count: bool, stack: bool,
                html: bool, pdf: bool):
    paths = list()
    table = load(file)
    for graph in iter_graphs(table, fields, count, stack):
        if html:
            paths.append(graph.write_html())
        if pdf:
            paths.append(graph.write_pdf())
    return paths


class PosBarGraph(OneSampleRefGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    @classmethod
    def xattr(cls):
        return POS_FIELD

    @classmethod
    @abstractmethod
    def yattr(cls):
        return ""

    @classmethod
    @abstractmethod
    def from_table(cls, *args, **kwargs):
        """ Construct a graph from a TableLoader. """
        return cls(*args, **kwargs)


class FieldPosBarGraph(PosBarGraph, ABC):
    """ Bar graph of a table with multiple types of bit fields, with one
    bar for each position in the sequence. """

    def __init__(self, *args, field_codes: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.codes = field_codes

    @classmethod
    @abstractmethod
    def ycount(cls):
        """ Whether the y-axis represents counts. """
        return False

    @classmethod
    def yattr(cls):
        return "Count" if cls.ycount() else "Fraction"

    @classmethod
    def get_table_field(cls, table: CountTable, code: str):
        """ Load one field from the table. """
        return (table.get_field_count(code) if cls.ycount()
                else table.get_field_frac(code).round(PRECISION))

    def title(self):
        fs = '/'.join(sorted(CountTable.FIELD_CODES[c] for c in self.codes))
        return f"{self.yattr()}s of {fs} bases in {self.sample} over {self.ref}"

    def graph_filename(self):
        sort_codes = "".join(sorted(self.codes))
        by = "base" if self.data_type() is pd.Series else "field"
        return f"{self.yattr()}_{sort_codes}_by_{by}".lower()

    @classmethod
    @abstractmethod
    def get_table_data(cls, table: RelPosTableLoader | MaskPosTableLoader,
                       field_codes: str) -> pd.DataFrame | pd.Series:
        """ Load all fields given by `field_codes` from the table. """
        return

    @classmethod
    def from_table(cls, table: RelPosTableLoader | MaskPosTableLoader,
                   field_codes: str):
        return cls(data=cls.get_table_data(table, field_codes), seq=table.seq,
                   out_dir=table.out_dir, sample=table.sample, ref=table.ref,
                   field_codes=field_codes)


class CountFieldPosBarGraph(FieldPosBarGraph, ABC):
    """ FieldPosBarGraph where each bar represents a count of reads. """

    @classmethod
    def ycount(cls):
        return True


class FractionFieldPosBarGraph(FieldPosBarGraph, ABC):
    """ FieldPosBarGraph where each bar represents a fraction of reads. """

    @classmethod
    def ycount(cls):
        return False


class OnePosBarGraph(PosBarGraph, ABC):
    """ Bar graph wherein each bar represents a base in a sequence. """

    @classmethod
    def data_type(cls):
        return pd.Series

    @classmethod
    def cmap_type(cls):
        return SeqColorMap

    def traces(self):
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
                                     marker_color=self.colors[base],
                                     hovertext=hovertext,
                                     hoverinfo="text"))
        return traces


class StackPosBarGraph(PosBarGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    outcomes for a position in a sequence. """

    @classmethod
    def data_type(cls):
        return pd.DataFrame

    @classmethod
    def cmap_type(cls):
        return RelColorMap

    def traces(self):
        traces = list()
        # Construct a trace for each field.
        for field, vals in self.data.items():
            # Define the text shown on hovering over a bar.
            hovertext = [f"{chr(base)}{x} {field}: {y}"
                         for base, (x, y) in zip(self.seq, vals.items(),
                                                 strict=True)]
            # Create a trace comprising all bars for this field.
            traces.append(go.Bar(name=field, x=vals.index, y=vals,
                                 marker_color=self.colors[field],
                                 hovertext=hovertext,
                                 hoverinfo="text"))
        return traces

    @cached_property
    def figure(self):
        fig = super().figure
        # Stack the bars at each position.
        fig.update_layout(barmode="stack")
        return fig


class SubplotPosBarGraph(PosBarGraph, ABC):

    @classmethod
    def data_type(cls):
        return pd.DataFrame

    @classmethod
    def cmap_type(cls):
        return SeqColorMap

    @property
    @abstractmethod
    def data_cols(self):
        """ Columns of the data to use. """
        return pd.Index()

    def traces(self):
        traces = dict()
        # Costruct a list of traces for each subplot.
        for col in self.data_cols:
            traces[col] = list()
            # Construct a trace for each type of base.
            for base in BASES:
                # Find the non-missing value at every base of that type.
                vals = self.data.loc[self.seqarr == base, col].dropna()
                # Check if there are any values to graph.
                if vals.size > 0:
                    # Define the text shown on hovering over a bar.
                    hovertext = [f"{chr(base)}{x}: {y}"
                                 for x, y in vals.items()]
                    # Create a trace comprising all bars for this base type.
                    traces[col].append(go.Bar(name=chr(base),
                                              x=vals.index, y=vals,
                                              marker_color=self.colors[base],
                                              hovertext=hovertext,
                                              hoverinfo="text"))
        return traces

    def _figure_init(self):
        return make_subplots(rows=self.data_cols.size, cols=1,
                             subplot_titles=self.data_cols)

    @cached_property
    def figure(self):
        fig = super().figure
        for row, (col, traces) in enumerate(self.traces().items(), start=1):
            for trace in traces:
                fig.add_trace(trace, row=row, col=1)
        return fig


class OneFieldPosBarGraph(FieldPosBarGraph, OnePosBarGraph, ABC):
    """ FieldPosBarGraph with a single bar for each sequence position. """

    @classmethod
    def get_table_data(cls, table: RelPosTableLoader | MaskPosTableLoader,
                       field_codes: str) -> pd.Series:
        if len(field_codes) != 1:
            raise ValueError(
                f"Expected 1 field code, but got {len(field_codes)}")
        return cls.get_table_field(table, field_codes[0])


class StackFieldPosBarGraph(FieldPosBarGraph, StackPosBarGraph, ABC):
    """ FieldPosBarGraph with stacked bars for each sequence position. """

    @classmethod
    def get_table_data(cls, table: RelPosTableLoader | MaskPosTableLoader,
                       field_codes: str):
        data = dict(((series := cls.get_table_field(table, code)).name, series)
                    for code in field_codes)
        return pd.DataFrame.from_dict(data)


class FullPosBarGraph(PosBarGraph, ABC):
    """ Bar graph of the positions in the full reference sequence. """

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.GraphSeg,)


class SectPosBarGraph(OneSampleSectGraph, PosBarGraph, ABC):
    """ Bar graph of the positions in a section. """

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.GraphSeg,)


class MaskPosBarGraph(FieldPosBarGraph, SectPosBarGraph, ABC):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """

    def title(self):
        return f"{super().title()}, section {self.sect}"

    @classmethod
    def from_table(cls, table: MaskPosTableLoader, field_codes: str):
        return cls(data=cls.get_table_data(table, field_codes), seq=table.seq,
                   out_dir=table.out_dir, sample=table.sample, ref=table.ref,
                   sect=table.sect, field_codes=field_codes)


class OneFractionRelatePosBarGraph(FractionFieldPosBarGraph,
                                   FullPosBarGraph,
                                   OneFieldPosBarGraph):
    """ """


class OneCountRelatePosBarGraph(CountFieldPosBarGraph,
                                FullPosBarGraph,
                                OneFieldPosBarGraph):
    """ """


class StackFractionRelatePosBarGraph(FractionFieldPosBarGraph,
                                     FullPosBarGraph,
                                     StackFieldPosBarGraph):
    """ """


class StackCountRelatePosBarGraph(CountFieldPosBarGraph,
                                  FullPosBarGraph,
                                  StackFieldPosBarGraph):
    """ """


class OneFractionMaskPosBarGraph(MaskPosBarGraph,
                                 FractionFieldPosBarGraph,
                                 OneFieldPosBarGraph):
    """ """


class OneCountMaskPosBarGraph(MaskPosBarGraph,
                              CountFieldPosBarGraph,
                              OneFieldPosBarGraph):
    """ """


class StackFractionMaskPosBarGraph(MaskPosBarGraph,
                                   FractionFieldPosBarGraph,
                                   StackFieldPosBarGraph):
    """ """


class StackCountMaskPosBarGraph(MaskPosBarGraph,
                                CountFieldPosBarGraph,
                                StackFieldPosBarGraph):
    """ """


class ClustPosBarGraph(SubplotPosBarGraph, SectPosBarGraph):
    """ Bar graph of a table of per-cluster mutation rates. """

    def __init__(self, *args, ks: Iterable[int], **kwargs):
        super().__init__(*args, **kwargs)
        self._ks = list(ks)

    @classmethod
    def yattr(cls):
        return "Mutation Rate"

    @property
    def sect(self):
        return self._sect

    def title(self):
        return (f"{self.yattr()}s per cluster in {self.sample} "
                f"over {self.ref}, section {self.sect}")

    @cached_property
    def data_cols(self):
        if not self._ks:
            # Return all columns.
            return self.data.columns.drop(SEQ_FIELD)
        # Return only the selected columns.
        return format_names_ks(self._ks)

    @property
    def data_labels(self):
        return parse_names(self.data_cols)

    def graph_filename(self):
        return f"clusters_{'-'.join(f'k{k}' for k in self._ks)}".lower()

    @classmethod
    def from_table(cls, table: ClusterPosTableLoader, ks: Iterable[int] = ()):
        return cls(data=table.data, seq=table.seq, out_dir=table.out_dir,
                   sample=table.sample, ref=table.ref, sect=table.sect, ks=ks)
