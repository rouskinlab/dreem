from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, Optional


OUTPUT_DIR = "output"
TEMP_DIR = "temp"

BASE = "base"
PARTITION = "partition"
MODULE = "module"

PARTITIONS = [OUTPUT_DIR, TEMP_DIR]
MODULES = ["alignment", "vectoring"]


class DreemPath(object):
    def __init__(self, name: str, parent: Optional[DreemPath] = None) -> None:
        self._name = Path(name)
        self._parent = parent
    
    @property
    def name(self):
        return self._name
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def path(self) -> Path:
        if self.parent is None:
            return self.name
        return self.parent.path.joinpath(self.name)


class BasePath(DreemPath):
    def __init__(self, path: str) -> None:
        super().__init__(path)


class PartitionPath(DreemPath):
    def __init__(self, name: str, parent: BasePath) -> None:
        assert name in PARTITIONS
        super().__init__(name, parent)


class ModulePath():
    def __init__(self, name: str, parent: PartitionPath) -> None:
        assert name in MODULES
        super().__init__(name, parent)



