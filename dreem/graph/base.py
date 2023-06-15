from abc import ABC, abstractmethod
from logging import getLogger
from functools import cached_property
from pathlib import Path

import pandas as pd
from plotly import graph_objects as go

from .color import get_cmap, ColorMap
from ..core import path
from ..core.sect import seq_to_int_array
from ..core.seq import DNA
from ..table.load import load, TableLoader

logger = getLogger(__name__)


def find_tables(tables: tuple[str, ...]):
    """ Return a file for each given file/directory of a table. """
    return path.find_files_multi(map(Path, tables), [path.MutTabSeg])


class GraphBase(ABC):
    def __init__(self, *, out_dir: Path, data: pd.DataFrame | pd.Series,
                 colors: str | None = None):
        self.out_dir = out_dir
        if not isinstance(data, self.data_type()):
            raise TypeError(f"{self.__class__.__name__} requires data of type "
                            f"'{self.data_type().__name__}', but got "
                            f"'{type(data).__name__}'")
        self.data = data
        self.colors = get_cmap(self.cmap_type(), colors)

    @abstractmethod
    def title(self):
        """ Title of the graph. """
        return ""

    @abstractmethod
    def xattr(self):
        """ Name of the attribute on the x-axis. """
        return ""

    @abstractmethod
    def yattr(self):
        """ Name of the attribute on the y-axis. """
        return ""

    @classmethod
    @abstractmethod
    def data_type(cls) -> type | tuple[type, ...]:
        """ Type of the data: either DataFrame, Series, or both. """
        return pd.DataFrame, pd.Series

    @classmethod
    @abstractmethod
    def cmap_type(cls):
        """ Type of the color map. """
        return ColorMap

    @abstractmethod
    def traces(self):
        """ Data traces. """
        return list()

    @abstractmethod
    def graph_filename(self):
        return ""

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Path segments. """
        return path.ModSeg,

    @abstractmethod
    def path_fields(self):
        """ Path fields. """
        return {path.TOP: self.out_dir, path.MOD: path.MOD_GRAPH,
                path.GRAPH: self.graph_filename()}

    def path(self, ext: str):
        """ Path to the output file of the graph. """
        return path.buildpar(*self.path_segs(), **self.path_fields(), ext=ext)

    def _figure_init(self):
        """ Initialize the figure. """
        return go.Figure(data=self.traces())

    def _figure_layout(self, fig: go.Figure):
        """ Update the figure's layout. """
        fig.update_layout(title=self.title(),
                          xaxis=dict(title=self.xattr()),
                          yaxis=dict(title=self.yattr()),
                          plot_bgcolor="#ffffff",
                          paper_bgcolor="#ffffff")
        fig.update_xaxes(linewidth=1,
                         linecolor="#000000",
                         autorange=True)
        fig.update_yaxes(gridcolor="#d0d0d0",
                         linewidth=1,
                         linecolor="#000000",
                         autorange=True)
        return fig

    @cached_property
    def figure(self):
        """ Figure object. """
        return self._figure_layout(self._figure_init())

    def write_html(self):
        """ Write the graph to an HTML file. """
        file = self.path(ext=path.HTML_EXT)
        if file.is_file():
            logger.warning(f"File exists: {file}")
        else:
            self.figure.write_html(file)
        return file

    def write_pdf(self):
        """ Write the graph to a PDF file. """
        file = self.path(ext=path.PDF_EXT)
        if file.is_file():
            logger.warning(f"File exists: {file}")
        else:
            self.figure.write_image(file)
        return file


class OneSampleGraph(GraphBase, ABC):
    """ Graph of one sample. """

    @property
    @abstractmethod
    def sample(self):
        """ Name of the sample. """
        return ""

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Path segments. """
        return super().path_segs() + (path.SampSeg,)

    def path_fields(self):
        return {**super().path_fields(), path.SAMP: self.sample}


class OneRefGraph(GraphBase, ABC):
    """ Graph of one reference sequence.  """

    @property
    @abstractmethod
    def ref(self):
        """ Name of the reference sequence. """
        return ""

    @property
    @abstractmethod
    def seq(self):
        """ Reference sequence as a DNA object. """
        return DNA(b"")

    @cached_property
    def seqarr(self):
        """ Reference sequence as an array of 8-bit integers. """
        return seq_to_int_array(self.seq)


class OneSectGraph(OneRefGraph):
    """ Graph of one section of one reference sequence. """

    @property
    @abstractmethod
    def sect(self):
        """ Name of the section of the reference sequence. """
        return ""


class OneSampleRefGraph(OneSampleGraph, OneRefGraph, ABC):

    def __init__(self, *args, sample: str, ref: str, seq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        if len(seq) != self.data.index.size:
            raise ValueError(f"Got different lengths of seq ({len(seq)}) and "
                             f"data ({self.data.index.size})")
        self._sample = sample
        self._ref = ref
        self._seq = seq

    @property
    def sample(self):
        return self._sample

    @property
    def ref(self):
        return self._ref

    @property
    def seq(self):
        return self._seq

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Path segments. """
        return super().path_segs() + (path.RefSeg,)

    def path_fields(self):
        return {**super().path_fields(), path.REF: self.ref}


class OneSampleSectGraph(OneSampleRefGraph, OneSectGraph, ABC):

    def __init__(self, *args, sect: str, **kwargs):
        super().__init__(*args, **kwargs)
        self._sect = sect

    @property
    def sect(self):
        return self._sect

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Path segments. """
        return super().path_segs() + (path.SectSeg,)

    def path_fields(self):
        return {**super().path_fields(), path.SECT: self.sect}


class GraphWriter(ABC):
    """ Write the proper types of graphs for a given table. """

    def __init__(self, table_file: Path):
        self.table_file = table_file

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return load(self.table_file)

    @abstractmethod
    def iter(self, fields: str, count: bool, stack: bool):
        """ Yield every graph for the table. """
        yield GraphBase(out_dir=Path(), data=pd.DataFrame())

    def write(self, fields: str, count: bool, stack: bool,
              html: bool, pdf: bool):
        """ Generate and write every graph for the table. """
        # Get the path for every graph.
        paths = list()
        for graph in self.iter(fields, count, stack):
            if html:
                paths.append(graph.write_html())
            if pdf:
                paths.append(graph.write_pdf())
        return paths
