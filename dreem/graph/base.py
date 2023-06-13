from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path

import pandas as pd
from plotly import graph_objects as go

from .color import get_cmap, ColorMap
from ..core import path
from ..table.load import load, PosTableLoader


def load_tables(tables: tuple[str, ...]):
    """ Return a TableLoader for each given path to a table. """
    return list(map(load, path.find_files_multi(map(Path, tables),
                                                [path.MutTabSeg])))


def load_pos_tables(tables: tuple[str, ...]):
    """ Return a PosTableLoader for each given path to a table that is
    indexed by position. """
    return list(filter(lambda table: isinstance(table, PosTableLoader),
                       load_tables(tables)))


class GraphBase(ABC):
    def __init__(self, *, out_dir: Path,
                 data: pd.DataFrame | pd.Series,
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

    @cached_property
    def figure(self):
        """ Figure object. """
        fig = go.Figure(data=self.traces())
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

    def write_html(self):
        file = self.path(ext=path.HTML_EXT)
        self.figure.write_html(file)
        return file

    def write_pdf(self):
        file = self.path(ext=path.PDF_EXT)
        self.figure.write_image(file)
        return file

    def write_png(self):
        file = self.path(ext=path.PNG_EXT)
        self.figure.write_image(file)
        return file
