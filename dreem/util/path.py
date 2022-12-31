from __future__ import annotations
import os
from pathlib import Path
import re
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
    def __init__(self, segment: str, require: bool=False) -> None:
        self._segment = segment
        if require:
            self.assert_exists()
    
    @property
    def path(self) -> Path:
        raise NotImplementedError
    
    @property
    def name(self) -> str:
        raise NotImplementedError
    
    @property
    def is_dir(self) -> bool:
        return True

    @property
    def is_file(self) -> bool:
        return not self.is_dir
    
    @property
    def exists(self) -> bool:
        if self.is_dir:
            return self.path.is_dir()
        if self.is_file:
            return self.path.is_file()
        assert False, f"Path '{self}' is neither directory nor file."
    
    def assert_exists(self):
        if not self.exists:
            raise FileNotFoundError(self.path)
    
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
        
        Arguments
        ---------
        path: str
            Path (absolute or relative) of the top-level directory into which
            all outputs from DREEM are written
        
        require: bool (default: False)
            Whether the path must exist at the time this function is called
        
        Returns
        -------
        None
        """
        super().__init__(sanitize(path), require)
    
    @property
    def path(self):
        return Path(self._segment)


class PathSubPart(PathPart):
    def __init__(self, parent: PathPart, name: str, require: bool=False):
        """
        Initialize instance of PathSubPart to represent directories and files
        within a top-level directory. This class represents directories, and
        the subclass PathFilePart represents files.

        Arguments
        ---------
        parent: PathLevel
            The immediate parent of the directory represented by this instance
        
        name: str
            Base name of the directory represented by this instance (i.e. the
            last component of the path, everything after the last separator).
            May not contain any path separators.
        
        require: bool (default: False)
            Whether the full path (including the final component, name) must
            exist at the time this instance is initialized
        
        Returns
        -------
        None

        Raises
        ------
        ValueError
            if name contains any path separators or parent does not represent
            the path of a directory (whether or not that directory exists)
        """
        if os.sep in name:
            raise ValueError(f"name ({name}) may not contain '{os.sep}'")
        if not parent.is_dir:
            raise ValueError(f"{parent} is not a directory: cannot be parent.")
        self._parent = parent
        super().__init__(name, require)
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def name(self):
        return self._segment
    
    @property
    def path(self):
        return self.parent.path.joinpath(self.name)

    def __getattribute__(self, name: str) -> Any:
        """
        Get an attribute of a PathSubPart instance. If the instance does not
        have the attribute, an AttributeError will be raised and caught.
        Then, the instance will check if its parent has the attribute.
        This search will continue up the chain of parents until either the
        attribute is found or the search reaches a TopDir, which neither
        implements this custom __getattribute__ method nor has a parent
        directory and will thus raise an uncaught AttributeError.
        This strategy of checking the parent for the attribute if the instance
        does not define it is designed to enable child directories to access
        attributes of their parents: for example, so a directory representing
        a region in a reference can determine which reference and sample the
        region belongs to.

        Warning: Using this method improperly can cause infinite recursion.
        
        Arguments
        ---------
        name: str
            Name of the attribute
        
        Returns
        -------
        any
            Value of the attribute for this instance, or if not found, for the
            closest parent with the attribute.
        
        Raises
        ------
        AttributeError
            if neither this instance nor any of its ancestors have the attribute
        """
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return self.parent.__getattribute__(name)


class DreemPath(object):
    parent_type: type = type(None)
    
    @classmethod
    def parse(cls, path: str, require: bool) -> DreemPath:
        raise NotImplementedError


class DreemTopDir(TopDir, DreemPath):
    def __init__(self, path: str, require: bool=False) -> None:
        super().__init__(path, require)

    @property
    def top(self):
        return self

    @classmethod
    def parse(cls, path: str, require: bool=True):
        return cls(path, require)


class DreemSubPath(PathSubPart, DreemPath):
    def __init__(self, parent: PathPart, name: str, require: bool=False):
        if not isinstance(parent, self.parent_type):
            raise ValueError(
                f"Parent of {type(self)} must be {self.parent_type}, "
                f"but got {type(parent)}.")
        super().__init__(parent, name, require)
    
    @classmethod
    def parse(cls, path: str, require: bool=True):
        """
        Parse a raw string and return a DreemPath.
        """
        head, tail = os.path.split(sanitize(path))
        parent = cls.parent_type.parse(head, require)
        return cls(parent, tail, require)


class DreemFile(DreemSubPath):
    @property
    def is_dir(self):
        return False


class DreemOutDir(DreemSubPath):
    parent_type = DreemTopDir

    def __init__(self, parent: parent_type, name: str, require: bool=False):
        if name not in (OUTPUT_DIR, TEMP_DIR):
            raise ValueError(f"Invalid output directory name: '{name}'")
        super().__init__(parent, name, require)
    
    @property
    def is_temp(self):
        return self.name == TEMP_DIR


class DreemModuleDir(DreemSubPath):
    parent_type = DreemOutDir

    def __init__(self, parent: parent_type, name: str, require: bool=False):
        if name not in MODULES:
            raise ValueError(f"Invalid module name: '{name}'")
        super().__init__(parent, name, require)
    
    @property
    def module(self):
        return self


class DreemSampleDir(DreemSubPath):
    parent_type = DreemModuleDir

    @property
    def sample(self):
        return self


class DreemRefDir(DreemSubPath):
    parent_type = DreemSampleDir

    @property
    def ref(self):
        return self


class DreemRegionDir(DreemSubPath):
    parent_type = DreemRefDir
    name_format = "{first}-{last}"
    name_pattern = re.compile("^([0-9]+)-([0-9]+)$")

    def __init__(self, parent: PathPart, name: str, require: bool=False):
        super().__init__(parent, name, require)
        self._first, self._last = self.extract_bounds(self.name_pattern, name)

    @staticmethod
    def extract_bounds(pattern: re.Pattern[str], name: str):
        if not (match := pattern.match(name)):
            raise ValueError(f"Invalid region name: '{name}'")
        first, last = map(int, match.groups())
        if not 1 <= first <= last:
            raise ValueError(f"Invalid region bounds: {first}-{last}")
        return first, last
    
    @classmethod
    def assemble(cls, parent: PathPart, first: int, last: int):
        name = cls.name_format.format(first=first, last=last)
        return cls(parent, name)

    @property
    def region(self):
        """ Enable any children of this directory to access its attributes """
        return self
    
    @property
    def first(self):
        """ The first position in the region (1-indexed, inclusive) """
        return self._first
    
    @property
    def last(self):
        """ The last position in the region (1-indexed, inclusive) """
        return self._last
    
    def get_report_file(self):
        return DreemMpReportFile.assemble(self.parent, self.first, self.last)


class DreemMpReportFile(DreemFile):
    parent_type = DreemRefDir
    name_format = "{first}-{last}_report.txt"
    name_pattern = re.compile("^([0-9]+)-([0-9]+)_report.txt$")

    def __init__(self, parent: PathPart, name: str, require: bool=False):
        super().__init__(parent, name, require)
        if self.is_temp:
            raise ValueError(f"Report cannot be in temporary directory.")
        if (mod := self.module.name) != MOD_VEC:
            raise ValueError(f"Module must be '{MOD_VEC}', but got '{mod}'.")
        self._first, self._last = DreemRegionDir.extract_bounds(
            self.name_pattern, name)

    @classmethod
    def assemble(cls, parent: PathPart, first: int, last: int):
        name = cls.name_format.format(first=first, last=last)
        return cls(parent, name)
    
    @property
    def first(self):
        """ The first position in the region (1-indexed, inclusive) """
        return self._first

    @property
    def last(self):
        """ The last position in the region (1-indexed, inclusive) """
        return self._last
    
    def get_region_dir(self):
        return DreemRegionDir.assemble(self.parent, self.first, self.last)
