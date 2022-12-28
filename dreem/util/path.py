import os
from pathlib import Path, PurePath
from typing import Any, Union, Dict


# Modules
MOD_DMX = "demultiplexing"
MOD_ALN = "alignment"
MOD_VEC = "vectoring"
MOD_CLS = "clustering"
MOD_AGG = "aggregation"


OUTPUT_DIR = "output"
TEMP_DIR = "temp"


# Path functions

def switch_directory(old_path: str, new_dir: str):
    return os.path.join(new_dir, os.path.basename(old_path))


def try_remove(file: str):
    try:
        os.remove(file)
    except OSError:
        pass


class DreemPath(object):
    def __init__(self, path: Union[str, PurePath]) -> None:
        """
        Initialize a DreemPath, the base class of paths in DREEM.
        """
        super().__init__()
        self._path = Path(path)
        self._setup()
    
    def _setup(self):
        raise NotImplementedError
    
    @property
    def path(self):
        return self._path
    
    @property
    def is_file(self):
        return self.path.is_file()
    
    @property
    def is_dir(self):
        return self.path.is_dir()
    
    def mkdir(self, *args, **kwargs):
        return self.path.mkdir(*args, **kwargs)
    
    def __str__(self) -> str:
        return str(self.path)


class OutputDir(DreemPath):
    def __init__(self, out_dir: str) -> None:
        super().__init__(out_dir)
    
    def _setup(self):
        self.mkdir(parents=True, exist_ok=True)


class ChildDir(DreemPath):
    def __init__(self, parent: DreemPath, name: str) -> None:
        super().__init__(PurePath(parent.path).joinpath(name))
        self._parent = parent
        self._name = name
    
    def _setup(self):
        self.mkdir(parents=False, exist_ok=True)
    
    @property
    def parent(self) -> DreemPath:
        return self._parent
    
    @property
    def name(self) -> str:
        return self._name
    
    def __getattr__(self, name: str) -> Any:
        return self.parent.__getattr__(name)


# Module directories

class ModuleDir(ChildDir):
    def __init__(self, out_dir: OutputDir, module: str) -> None:
        super().__init__(out_dir, module)
        self._out_dir = out_dir
    
    @property
    def out_dir(self) -> OutputDir:
        return self._out_dir


class DemultiDir(ModuleDir):
    def __init__(self, out_dir: OutputDir) -> None:
        super().__init__(out_dir, MOD_DMX)


class AlignDir(ModuleDir):
    def __init__(self, out_dir: OutputDir) -> None:
        super().__init__(out_dir, MOD_ALN)


class VectorDir(ModuleDir):
    def __init__(self, out_dir: OutputDir) -> None:
        super().__init__(out_dir, MOD_VEC)


class ClusterDir(ModuleDir):
    def __init__(self, out_dir: OutputDir) -> None:
        super().__init__(out_dir, MOD_CLS)


class AggregateDir(ModuleDir):
    def __init__(self, out_dir: OutputDir) -> None:
        super().__init__(out_dir, MOD_AGG)


# Sample directories

class SampleDir(ChildDir):
    def __init__(self, mod_dir: ModuleDir, sample: str) -> None:
        super().__init__(mod_dir, sample)
        self._module = mod_dir
    
    @property
    def module(self):
        return self._module


# Reference directories

class SampRefDir(ChildDir):
    def __init__(self, samp_dir: SampleDir, ref: bytes) -> None:
        super().__init__(samp_dir, ref.decode())
        self._sample = samp_dir
    
    @property
    def sample(self) -> SampleDir:
        return self._sample

    @property
    def name(self) -> bytes:
        return self._name.encode()


# Region directories

class SampRefRegionDir(ChildDir):
    def __init__(self, ref_dir: SampRefDir, first: int, last: int) -> None:
        self._first = first
        self._last = last
        super().__init__(ref_dir, self.region)
    
    @property
    def region(self):
        return f"{self._first}-{self._last}"



MOD_DIRS: Dict[str, ModuleDir] = {
    MOD_DMX: DemultiDir,
    MOD_ALN: AlignDir,
    MOD_VEC: VectorDir,
    MOD_CLS: ClusterDir,
    MOD_AGG: AggregateDir
}


def new_dir(out: str, module: str = "", sample: str = "", ref: bytes = b"",
            region_first: int = 0, region_last: int = 0):
    current = OutputDir(out)
    if module:
        try:
            dir_cls = MOD_DIRS[module]
        except KeyError:
            raise ValueError(f"Invalid module name: '{module}'")
        else:
            current = dir_cls(current)
    if sample:
        if module:
            samp_dir = SampleDir(current, sample)
        else:
            raise ValueError("Cannot give sample without module.")



