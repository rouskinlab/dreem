from functools import cached_property
import os
from pathlib import Path, PurePath
from typing import Any, Union, Dict


# Modules
MOD_DMX = "demultiplexing"
MOD_ALN = "alignment"
MOD_VEC = "vectoring"
MOD_CLS = "clustering"
MOD_AGG = "aggregation"
MODULES = (MOD_DMX, MOD_ALN, MOD_VEC, MOD_CLS, MOD_AGG)


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

def switch_directory(old_path: str, new_dir: str):
    return os.path.join(new_dir, os.path.basename(old_path))


def try_remove(file: str):
    try:
        os.remove(file)
    except OSError:
        pass


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
    def has_first(self):
        if self.first == DEF_REGION_FIRST:
            return False
        if self.first <= 0:
            raise ValueError(f"first ({self.first}) must be positive integer.")
        return True
    
    @cached_property
    def has_last(self):
        if self.last == DEF_REGION_LAST:
            return False
        if self.last <= 0:
            raise ValueError(f"last ({self.last}) must be positive integer.")
        return True
    
    @cached_property
    def has_first_and_last(self):
        if self.has_first != self.has_last:
            raise ValueError("Must give both or neither of first/last.")
        return self.has_first

    @cached_property
    def has_region_dir(self):
        if (has_fl := self.has_first_and_last) and not self.has_ref_dir:
            raise ValueError("Got first and last without ref_name.")
        return has_fl and self.mv_report == DEF_MV_REPORT
    
    @property
    def is_mv_report_file(self):
        if self.mv_report == DEF_MV_REPORT:
            return False
        if not self.has_first_and_last:
            raise ValueError("Got mv_report without first and last.")
        return True
    
    @property
    def is_mv_batch_file(self):
        if self.mv_batch == DEF_MV_BATCH:
            return False
        if not self.has_region_dir:
            raise ValueError("Got mv_batch without first and last.")
        if self.mv_batch < 0:
            raise ValueError(f"mv_batch ({self.mv_batch}) must be "
                             "a non-negative integer.")
        return True


    def makedirs(self):
        if self.of_dir:
            os.makedirs(self.path)
        else:
            raise ValueError(f"Not a directory: '{self}'")
    
    @property
    def are_files(self):
        return {
            "mv_report": self.is_mv_report_file,
            "mv_batch": self.is_mv_batch_file
        }
    
    @property
    def of_file(self):
        file_sum = sum(self.are_files.values())
        if file_sum > 1:
            raise ValueError("Path specified as more than one file type.")
        return bool(file_sum)
    
    @property
    def of_dir(self):
        return not self.of_file
    
    @property
    def exists(self):
        return self.path.exists()
    
    def validate(self):
        _ = self.path
        _ = self.of_file
    
    def __str__(self) -> str:
        return str(self.path)
