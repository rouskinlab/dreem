import os
import shutil

from dreem.util.path import TEMP_DIR


# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
FASTQ_REC_LENGTH = 4
DEFAULT_MIN_MAPQ = 30



class NgsFileBase(object):
    _operation_dir = ""

    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str) -> None:
        self._root_dir = root_dir
        self._ref_file = ref_file
        self._sample = sample
    
    @property
    def operation_dir(self):
        if not self._operation_dir:
            raise NotImplementedError
        return self._operation_dir
    
    @property
    def root_dir(self):
        return self._root_dir
    
    @property
    def ref_file(self):
        return self._ref_file
    
    @property
    def ref_prefix(self):
        return os.path.splitext(self.ref_file)[0]
    
    @property
    def ref_filename(self):
        return os.path.basename(self.ref_prefix)
    
    @property
    def sample(self):
        return self._sample
    
    def _get_dir(self, dest: str):
        return os.path.join(self.root_dir, dest,
                            self.operation_dir, self.sample)
    
    @property
    def output_dir(self):
        return self._get_dir(TEMP_DIR)
    
    def _make_output_dir(self):
        os.makedirs(self.output_dir, exist_ok=True)
    
    def run(self, *args, **kwargs):
        raise NotImplementedError
    
    def clean(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)
