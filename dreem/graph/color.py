from abc import ABC, abstractmethod
from functools import cache
from inspect import getmembers
from sys import modules
from typing import Any, Hashable

from ..core.seq import A_INT, C_INT, G_INT, T_INT
from ..table.base import CountTable


class ColorMap(ABC):
    """ Color map for a graph. """

    def __init__(self, name: str, **kwargs):
        self.name = name
        self._colors = self._set_colors(**kwargs)

    @abstractmethod
    def _set_colors(self, **kwargs):
        return dict(**kwargs)

    def get(self, item: Hashable, default: Any | None = None):
        return self._colors.get(item, default)

    def __getitem__(self, item: Hashable):
        return self._colors[item]


class SeqColorMap(ColorMap):
    """ Color map for bases A, C, G, and T. """

    def __init__(self, name: str, a: str, c: str, g: str, t: str):
        super().__init__(name, a=a, c=c, g=g, t=t)

    def _set_colors(self, *, a: str, c: str, g: str, t: str):
        return {A_INT: a, C_INT: c, G_INT: g, T_INT: t}


class RelColorMap(ColorMap):
    """ Color map for relationships. """

    def __init__(self, name: str, b: str, n: str, r: str, m: str,
                 d: str, i: str, s: str, a: str, c: str, g: str, t: str):
        super().__init__(name, b=b, n=n, r=r, m=m,
                         d=d, i=i, s=s, a=a, c=c, g=g, t=t)

    def _set_colors(self, **kwargs):
        colors = dict()
        for key, field in CountTable.FIELD_CODES.items():
            colors[field] = kwargs.pop(key)
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {kwargs}")
        return colors


basic = SeqColorMap("basic", a="#FF0000", c="#0000FF", g="#FFC000", t="#008000")
water = SeqColorMap("water", a="#A15252", c="#3D427D", g="#E3CC7B", t="#76B887")
earth = SeqColorMap("earth", a="#D17777", c="#464EA6", g="#E3CC7B", t="#336140")
steel = SeqColorMap("steel", a="#663328", c="#716B80", g="#91B8AC", t="#D9D5B4")

crayons = RelColorMap("crayons", b="#531B93", n="#009193", r="#73FCD6",
                      m="#FF40FF", d="#FF9300", i="#00F900", s="#7A81FF",
                      a="#941100", c="#0096FF", g="#FFFC79", t="#4F8F00")


DEFAULTS: dict[type[ColorMap], str] = {
    RelColorMap: "crayons",
    SeqColorMap: "earth",
}


@cache
def get_colormaps(cmap_class: type[ColorMap]):
    """ Return a dict of all color maps of a given class. """
    colormaps: dict[str, cmap_class] = dict()
    for _, cmap in getmembers(modules[__name__],
                              lambda item: isinstance(item, cmap_class)):
        if cmap.name in colormaps:
            raise ValueError(f"Duplicate {cmap_class.__name__}: '{cmap.name}'")
        colormaps[cmap.name] = cmap
    if (default := DEFAULTS[cmap_class]) not in colormaps:
        raise ValueError(
            f"Default {cmap_class.__name__} '{default}' does not exist")
    return colormaps


def get_cmap(cmap_class: type[ColorMap], name: str | None = None):
    """ Get a color map of a given class by its name. """
    cmaps = get_colormaps(cmap_class)
    if not cmaps:
        raise ValueError(f"No color maps of class {cmap_class.__name__}")
    if name is None:
        # Use the default color map for the class.
        name = DEFAULTS[cmap_class]
    return cmaps[name]
