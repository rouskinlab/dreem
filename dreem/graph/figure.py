from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
from plotly import graph_objects as go

from ..core import path
from ..core.seq import DNA, BASES, A_INT, C_INT, G_INT, T_INT, seq_to_int_array
from ..table.load import POS_FIELD, PosTableLoader


# Define color encodings for each type of base.

class BaseColors(object):
    def __init__(self, a: str, c: str, g: str, t: str, default: str):
        self._colors = {A_INT: a, C_INT: c, G_INT: g, T_INT: t}
        self._default = default

    def get(self, item):
        return self._colors.get(item, self._default)

    def __getitem__(self, item):
        return self._colors[item]


base_colors = {
    "basic": BaseColors(a="#ff0000", c="#0000ff", g="#ffc000", t="#008000",
                        default="#808080"),
    "mellow": BaseColors(a="#C23C3C", c="#6478F5", g="#B5A73E", t="#42B38F",
                         default="bdbdbd")
}

PRECISION = 6  # number of decimal places for hover text


class GraphBase(ABC):
    def __init__(self, data: pd.DataFrame | pd.Series,
                 sample: str, ref: str):
        self.data = data
        self.sample = sample
        self.ref = ref

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

    @abstractmethod
    def traces(self):
        """ Data traces. """
        return list()

    @abstractmethod
    def path_segs(self):
        """ Path segments. """
        return path.ModSeg, path.SampSeg, path.RefSeg

    @abstractmethod
    def path_fields(self, out_dir: Path, ext: str):
        """ Path fields. """
        return {path.TOP: out_dir, path.MOD: path.MOD_GRAPH,
                path.SAMP: self.sample, path.REF: self.ref, path.EXT: ext}

    def path(self, out_dir: Path, ext: str):
        """ Path to the output file of the graph. """
        return path.build(*self.path_segs(), **self.path_fields(out_dir, ext))

    def figure(self):
        """ Figure object. """
        fig = go.Figure(data=self.traces())
        fig.update_layout(title=self.title,
                          xaxis=dict(title=self.xattr()),
                          yaxis=dict(title=self.yattr()),
                          plot_bgcolor="#ffffff",
                          paper_bgcolor="#ffffff")
        fig.update_xaxes(
            linewidth=1,
            linecolor="#000000",
            mirror=True,
            autorange=True
        )
        fig.update_yaxes(
            gridcolor="#d0d0d0",
            linewidth=1,
            linecolor="#000000",
            mirror=True,
            autorange=True
        )
        return fig


class SeqBarGraphBase(GraphBase, ABC):
    def __init__(self, title: str, data: pd.DataFrame | pd.Series, seq: DNA, *,
                 colors: str):
        super().__init__(title, data)
        self.seq = seq
        if len(self.seq) != data.index.size:
            raise ValueError(f"Got {len(self.seq)} bases in seq "
                             f"but {data.index.size} positions in data")
        self.colors = base_colors[colors]

    @classmethod
    def xattr(cls):
        return POS_FIELD

    @cached_property
    def seqarr(self):
        return seq_to_int_array(self.seq)


class SeqOneBarGraph(SeqBarGraphBase):
    def __init__(self, title: str, data: pd.Series, seq: DNA):
        super().__init__(title, data, seq, colors="mellow")

    def yattr(self):
        return self.data.name

    def traces(self):
        traces = list()
        # Construct a trace for each type of base.
        for base in BASES:
            # Find the positions at which the base is this type of base.
            base_pos = self.data.index[self.seqarr == base]
            if base_pos.size > 0:
                # Define the text shown on hovering over a bar.
                hovertext = [
                    f"{chr(base)}{pos}: {round(self.data.loc[pos], PRECISION)}"
                    for pos in base_pos
                ]
                # Create a trace comprising all bars for this base type.
                traces.append(go.Bar(
                    name=chr(base),
                    x=base_pos,
                    y=self.data.loc[base_pos],
                    marker_color=self.colors[base],
                    hovertext=hovertext,
                    hoverinfo="text",
                ))
        return traces

    @classmethod
    def from_table(cls, loader: PosTableLoader, yattr: str):
        return cls(yattr, loader.data[yattr], loader.seq)



def mutation_fraction(df, show_ci: bool = True) -> dict:
    assert len(df) == 1, "df must have only one row"
    mh = df.iloc[0]
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map

    traces, layouts = [], []
    mh['index_selected'] = [i + 1 for i in
                            range(len(mh['sequence']))]  # TODO[i + 1 for i in mh.index_selected] # index starts at 1
    mh_unrolled = pd.DataFrame(
        {'mut_rate': list(mh.sub_rate), 'base': list(mh.sequence), 'index_reset': list(range(len(mh.index_selected))),
         'index_selected': mh.index_selected, 'paired': list(mh.structure)})

    for bt in set(mh['sequence']):
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate': list(df_loc.mut_rate),
                                   'base': list(df_loc.base),
                                   'index': df_loc['index_selected'],
                                   'paired': [{'.': False, '(': True, ')': True}[s] for s in df_loc.paired]})
        traces.append(go.Bar(
            x=np.array(df_loc['index_reset']),
            y=np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text=hover_attr,
            hovertemplate=''.join(["<b>" + ha + ": %{text[" + str(i) + "]}<br>" for i, ha in enumerate(hover_attr)]),
        ))
        if show_ci:
            traces[-1].update(
                error_y=dict(
                    type='data',
                    symmetric=False,
                ))

    fig = go.Figure(data=traces)

    fig.update_layout(
        title=f"{mh['sample']} - {mh['reference']} - {mh['section']} - {mh['cluster']} - {mh['num_aligned']} reads",
        xaxis=dict(title="Position"),
        yaxis=dict(title="Mutation fraction", range=[0, 0.1]))

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )

    fig.update_xaxes(
        tickvals=mh_unrolled['index_reset'],
        ticktext=["%s %s" % ({'.': '(U)', '(': '(P)', ')': '(P)'}[x], str(y)) for (x, y) in
                  zip(mh['structure'], mh['index_selected'])],
        tickangle=90,
        autorange=True
    )

    # make the background white
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': mh}
