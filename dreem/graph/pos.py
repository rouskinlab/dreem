from abc import ABC, abstractmethod
from functools import cached_property
import os
from pathlib import Path

from click import command
import pandas as pd
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from .base import GraphBase, load_pos_tables
from .color import RelColorMap, SeqColorMap
from ..core import docdef, path
from ..core.cli import opt_table, opt_max_procs, opt_parallel
from ..core.parallel import dispatch, as_list_of_tuples
from ..core.seq import BASES, DNA, seq_to_int_array
from ..table.base import CountTable
from ..table.load import (POS_FIELD, PosTableLoader, MaskPosTableLoader,
                          RelPosTableLoader, ClusterPosTableLoader)

# Number of digits to which to round decimals.
PRECISION = 6

params = [
    opt_table,
    opt_max_procs,
    opt_parallel,
]


@command(__name__.split(os.path.extsep)[-1], params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positional attributes. """
    return run(*args, **kwargs)


@docdef.auto()
def run(table: tuple[str, ...],
        yaxis: tuple[str, ...] = (), *,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    return dispatch(graph_table, max_procs, parallel, pass_n_procs=False,
                    args=as_list_of_tuples(load_pos_tables(table)))


def graph_table(table: PosTableLoader):
    if isinstance(table, ClusterPosTableLoader):
        graph = ClustPosGraph.from_table(table)
    else:
        return
        graph = StackFractionMaskFieldPosGraph.from_table(tab, "acgtdi")
    graph.write_html()


class PosGraphBase(GraphBase, ABC):
    """ Base class for bar graphs wherein each bar represents a position
    in a sequence. """

    def __init__(self, *, seq: DNA, sample: str, ref: str, **kwargs):
        super().__init__(**kwargs)
        if len(seq) != self.data.index.size:
            raise ValueError(f"Got different lengths of seq ({len(seq)}) and "
                             f"data ({self.data.index.size})")
        self.seq = seq
        self.sample = sample
        self.ref = ref

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Path segments. """
        return super().path_segs() + (path.SampSeg, path.RefSeg)

    @abstractmethod
    def path_fields(self):
        return {**super().path_fields(),
                path.SAMP: self.sample, path.REF: self.ref}

    @classmethod
    def xattr(cls):
        return POS_FIELD

    @cached_property
    def seqarr(self):
        return seq_to_int_array(self.seq)

    @classmethod
    @abstractmethod
    def yattr(cls):
        return ""

    @classmethod
    @abstractmethod
    def from_table(cls, *args, **kwargs):
        """ Construct a graph from a TableLoader. """
        return cls(*args, **kwargs)


class FieldPosGraph(PosGraphBase, ABC):
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


class CountFieldPosGraph(FieldPosGraph, ABC):
    """ FieldPosGraph where each bar represents a count of reads. """

    @classmethod
    def ycount(cls):
        return True


class FractionFieldPosGraph(FieldPosGraph, ABC):
    """ FieldPosGraph where each bar represents a fraction of reads. """

    @classmethod
    def ycount(cls):
        return False


class SingleBarPosGraph(PosGraphBase, ABC):
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


class StackPosGraph(PosGraphBase, ABC):
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


class SubplotPosGraph(PosGraphBase, ABC):

    @classmethod
    def data_type(cls):
        return pd.DataFrame

    @classmethod
    def cmap_type(cls):
        return SeqColorMap

    def traces(self):
        traces = dict()
        # Costruct a list of traces for each subplot.
        for col, series in self.data.items():
            traces[col] = list()
            # Construct a trace for each type of base.
            for base in BASES:
                # Find the non-missing value at every base of that type.
                vals = self.data.loc[self.seqarr == base].dropna()
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

    @cached_property
    def figure(self):
        fig = make_subplots(rows=self.data.columns.size, cols=1)
        for row, (col, traces) in enumerate(self.traces().items(), start=1):
            for trace in traces:
                fig.add_trace(trace)
        return fig


class SingleBarFieldPosGraph(FieldPosGraph, SingleBarPosGraph, ABC):
    """ FieldPosGraph with a single bar for each sequence position. """

    @classmethod
    def get_table_data(cls, table: RelPosTableLoader | MaskPosTableLoader,
                       field_codes: str) -> pd.Series:
        if len(field_codes) != 1:
            raise ValueError(
                f"Expected 1 field code, but got {len(field_codes)}")
        return cls.get_table_field(table, field_codes[0])


class StackFieldPosGraph(FieldPosGraph, StackPosGraph, ABC):
    """ FieldPosGraph with stacked bars for each sequence position. """

    @classmethod
    def get_table_data(cls, table: RelPosTableLoader | MaskPosTableLoader,
                       field_codes: str):
        data = dict(((series := cls.get_table_field(table, code)).name, series)
                    for code in field_codes)
        return pd.DataFrame.from_dict(data)


class SectPosGraph(PosGraphBase, ABC):
    """ Bar graph of the positions in a section. """

    @property
    @abstractmethod
    def sect(self):
        """ Name of the section. """
        return ""

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.SectSeg, path.GraphSeg)

    def path_fields(self):
        return {**super().path_fields(), path.SECT: self.sect}


class MaskFieldPosGraph(FieldPosGraph, SectPosGraph, ABC):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """

    def __init__(self, *args, sect: str, **kwargs):
        super().__init__(*args, **kwargs)
        self._sect = sect

    @property
    def sect(self):
        return self._sect

    def title(self):
        return f"{super().title()}, section {self.sect}"

    @classmethod
    def from_table(cls, table: MaskPosTableLoader, field_codes: str):
        return cls(data=cls.get_table_data(table, field_codes), seq=table.seq,
                   out_dir=table.out_dir, sample=table.sample, ref=table.ref,
                   sect=table.sect, field_codes=field_codes)


class SingleBarFractionMaskFieldPosGraph(MaskFieldPosGraph,
                                         FractionFieldPosGraph,
                                         SingleBarFieldPosGraph):
    """ """


class StackFractionMaskFieldPosGraph(MaskFieldPosGraph,
                                     FractionFieldPosGraph,
                                     StackFieldPosGraph):
    """ """


class ClustPosGraph(SubplotPosGraph, SectPosGraph):
    """ Bar graph of a table of per-cluster mutation rates. """

    def __init__(self, *args, sect: str, **kwargs):
        super().__init__(*args, **kwargs)
        self._sect = sect

    @classmethod
    def yattr(cls):
        return "Mutation Rate"

    @property
    def sect(self):
        return self._sect

    def title(self):
        return (f"{self.yattr()}s of bases in {self.sample} over {self.ref}, "
                f"section {self.sect}")

    def graph_filename(self):
        return f"clusters".lower()

    @classmethod
    def from_table(cls, table: ClusterPosTableLoader):
        return cls(data=table.data, seq=table.seq, out_dir=table.out_dir,
                   sample=table.sample, ref=table.ref, sect=table.sect)
