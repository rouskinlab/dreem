from __future__ import annotations
from functools import cached_property
import os
from pathlib import Path
from typing import Any, Optional, Callable


# Define names of module directories.
MOD_DMX = "demultiplexing"
MOD_ALN = "alignment"
MOD_VEC = "vectoring"
MOD_CLS = "clustering"
MOD_AGG = "aggregation"
MODULES = (MOD_DMX, MOD_ALN, MOD_VEC, MOD_CLS, MOD_AGG)


# Define default values for path components.
DEF_OUT_ROOT = os.getcwd()
DEF_MODULE = ""
DEF_SAMPLE = ""
DEF_REF_NAME = ""
DEF_REGION_FIRST = -1
DEF_REGION_LAST = -1
DEF_MV_BATCH = -1
DEF_MV_REPORT = False


OUTPUT_DIR = "output"
TEMP_DIR = "temp"


# Path functions

def sanitize(path: str):
    return os.path.realpath(os.path.normpath(os.path.abspath(path)))


def switch_directory(old_path: str, new_dir: str):
    return os.path.join(new_dir, os.path.basename(old_path))


def try_remove(file: str):
    try:
        os.remove(file)
    except OSError:
        pass


class PathPart(object):
    def __init__(self, node: str) -> None:
        self._node = node
    
    @property
    def path(self) -> Path:
        raise NotImplementedError
    
    @property
    def of_dir(self) -> bool:
        return True

    @property
    def of_file(self) -> bool:
        return not self.of_dir
    
    @property
    def exists(self) -> bool:
        if self.of_dir:
            return self.path.is_dir()
        if self.of_file:
            return self.path.is_file()
        assert False, f"Path '{self}' is neither directory nor file."
    
    def assert_exists(self):
        if not self.exists:
            raise OSError(f"Path '{self}' does not exist.")
    
    def __str__(self) -> str:
        return str(self.path)


class TopDir(PathPart):
    def __init__(self, path: str, require: bool = False) -> None:
        """
        Initialize instance of PathTopLevel to represent the top-level
        directory of a DREEM path. Top-level means the directory at and above
        which the file structure is not relevant to the operation of DREEM.
        For example, if your project's files are in the directory
          /home/me/projects/my_fav_rna
        and you tell DREEM to output all results to the directory
          /home/me/projects/my_fav_rna/results
        then the top-level directory would be
          /home/me/projects/my_fav_rna/results
        which would contain any of the following sub-directories
          /home/me/projects/my_fav_rna/results/demultiplexing
          /home/me/projects/my_fav_rna/results/alignment
          /home/me/projects/my_fav_rna/results/vectoring
          /home/me/projects/my_fav_rna/results/clustering
          /home/me/projects/my_fav_rna/results/aggregation
        
        * Arguments *
        path (str) -----> path (absolute or relative) of the top-level directory
                          into which all outputs from DREEM are written
        require (bool) -> whether the path must exist at the time this function
                          is called (default: False)
        
        * Returns *
        None
        """
        super().__init__(sanitize(path))
        if require:
            self.assert_exists()
    
    @property
    def path(self):
        return Path(self._node)


class PathSubPart(PathPart):
    def __init__(self, name: str, parent: PathPart, require: bool = False,
                 filter: Optional[Callable[[str, PathPart], None]] = None):
        """
        Initialize instance of PathSubPart to represent directories and files
        within a top-level directory. This class represents directories, and
        the subclass PathFilePart represents files.

        * Arguments *
        name (str) -------------> base name of the directory represented by this
                                  instance (i.e. the last component of the path,
                                  everything after the last slash). May contain
                                  path separator characters (and thus represent
                                  multiple directories), but this is discouraged
                                  as it may cause confusion.
        parent (PathLevel) -----> directory in which the base of the path is
                                  located: the immediate parent of the base path
        require (bool) ---------> whether the full path (including the final
                                  component, name) must exist at the time this
                                  instance is initialized (default: False).
        filter (function|None) -> (optional) validate the values of name and
                                  parent with the function passed to filter.
                                  If given, this function must accept name and
                                  path as arguments, return None if successful,
                                  and raise an exception if the validation
                                  fails (default: None = no validation).

        * Returns *
        None
        """
        super().__init__(name)
        if not parent.of_dir:
            raise ValueError(f"{parent} is not a directory: cannot be parent.")
        if filter:
            filter(name, parent)
        self._parent = parent
        if require:
            self.assert_exists()
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def name(self):
        return self._node
    
    @property
    def path(self):
        return self.parent.path.joinpath(self.name)


class DreemPath(object):
    parent_type: type = type(None)
    
    @classmethod
    def parse(cls, path: str) -> DreemPath:
        raise NotImplementedError


class DreemTopDir(TopDir, DreemPath):
    def __init__(self, path: str, require: bool = False) -> None:
        super().__init__(path, require)

    @property
    def top(self):
        return self

    @classmethod
    def parse(cls, path: str, require: bool = False):
        return cls(path, require)


class DreemSubPath(PathSubPart, DreemPath):
    def __init__(self, name: str, parent: PathPart, require: bool = False,
                 filter: Optional[Callable[[str, PathPart], None]] = None):
        if not isinstance(parent, self.parent_type):
            raise ValueError(
                f"Parent of {type(self)} must be {self.parent_type}, "
                f"but got {type(parent)}.")
        super().__init__(name, parent, require, filter)
    
    def __getattr__(self, attr: str) -> Any:
        return self.parent.__getattribute__(attr)
    
    @classmethod
    def parse(cls, path: str, require: bool = False):
        """
        Parse a raw string and return a DreemPath.
        """
        sanitized = sanitize(path)
        head, tail = os.path.split(sanitized)
        parent = cls.parent_type.parse(head, require)
        return cls(tail, parent, require=require)


class DreemFile(DreemPath):
    @property
    def of_dir(self):
        return False


class DreemModuleDir(DreemSubPath):
    parent_type = DreemTopDir

    def __init__(self, name: str, parent: parent_type,
                 require: bool = False) -> None:
        def filter(name: str, parent: self.parent_type):
            if name not in MODULES:
                raise ValueError(f"Invalid module name: '{name}'")
        super().__init__(name, parent, require, filter)
    
    @property
    def module(self):
        return self


class DreemPath(object):
    def __init__(self,
                 out_root: str = DEF_OUT_ROOT,
                 module: str = DEF_MODULE,
                 sample: str = DEF_SAMPLE,
                 ref_name: str = DEF_REF_NAME,
                 first: int = DEF_REGION_FIRST,
                 last: int = DEF_REGION_LAST,
                 mv_batch: int = DEF_MV_BATCH,
                 mv_report: bool = DEF_MV_REPORT,
                 ) -> None:
        """
        Initialize a DreemPath to manipulate file paths in DREEM.
        """
        self.out_root = out_root
        self.module = module
        self.sample = sample
        self.ref_name = ref_name
        self.first = first
        self.last = last
        self.mv_batch = mv_batch
        self.mv_report = mv_report
        self.validate()

    @cached_property
    def path(self) -> Path:
        p = Path(self.out_root)
        if self.has_module_dir:
            p = p.joinpath(self.module)
        if self.has_sample_dir:
            if not self.has_module_dir:
                raise ValueError("Got sample without module.")
            p = p.joinpath(self.sample)
        if self.has_ref_dir:
            p = p.joinpath(self.ref_name)
        if self.has_region_dir:
            p = p.joinpath(f"{self.first}-{self.last}")
        elif self.is_mv_report_file:
            return p.joinpath(f"{self.first}-{self.last}_report.txt")
        if self.is_mv_batch_file:
            return p.joinpath(f"{self.mv_batch}.orc")
        return p

    @cached_property
    def has_module_dir(self):
        if self.module == DEF_MODULE:
            return False
        if self.module not in MODULES:
            raise ValueError(f"Invalid module name: '{self.module}'")
        return True
    
    @cached_property
    def has_sample_dir(self):
        if self.sample == DEF_SAMPLE:
            return False
        if not self.has_module_dir:
            raise ValueError("Got sample without module.")
        return True

    @cached_property
    def has_ref_dir(self):
        if self.ref_name == DEF_REF_NAME:
            return False
        if not self.has_sample_dir:
            raise ValueError("Got ref_name without sample.")
        return True
    
    @cached_property
    def _has_first(self):
        if self.first == DEF_REGION_FIRST:
            return False
        if self.first <= 0:
            raise ValueError(f"first ({self.first}) must be positive integer.")
        return True
    
    @cached_property
    def _has_last(self):
        if self.last == DEF_REGION_LAST:
            return False
        if self.last <= 0:
            raise ValueError(f"last ({self.last}) must be positive integer.")
        return True
    
    @cached_property
    def has_first_and_last(self):
        if self._has_first != self._has_last:
            raise ValueError("Must give both or neither of first/last.")
        if self._has_first and not self.has_ref_dir:
            raise ValueError("Got first and last without ref_name.")
        return self._has_first

    @cached_property
    def has_region_dir(self):
        return self.has_first_and_last and not self.is_mv_report_file
    
    @cached_property
    def is_mv_report_file(self):
        if self.mv_report == DEF_MV_REPORT:
            return False
        if not self.has_first_and_last:
            raise ValueError("Got mv_report without first and last.")
        return True
    
    @cached_property
    def is_mv_batch_file(self):
        if self.mv_batch == DEF_MV_BATCH:
            return False
        if not self.has_first_and_last:
            raise ValueError("Got mv_batch without first and last.")
        if self.mv_batch < 0:
            raise ValueError(f"mv_batch ({self.mv_batch}) must be "
                             "a non-negative integer.")
        return True

    def makedirs(self):
        if self.represents_dir:
            os.makedirs(self.path)
        else:
            raise ValueError(f"Not a directory: '{self}'")
    
    @property
    def is_each_file(self):
        return {
            "mv_report": self.is_mv_report_file,
            "mv_batch": self.is_mv_batch_file
        }
    
    @property
    def represents_file(self):
        file_sum = sum(self.is_each_file.values())
        if file_sum > 1:
            raise ValueError("Path specified as more than one file type.")
        return bool(file_sum)
    
    @property
    def represents_dir(self):
        return not self.represents_file
    
    @property
    def exists(self):
        return self.path.exists()
    
    def validate(self):
        _ = self.path
        _ = self.represents_dir
    
    def __str__(self) -> str:
        return str(self.path)


def parse_mv_report(cls, path: str):
    """
    Parse a given path and return a DreemPath instance.
    """