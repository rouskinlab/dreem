import os


# Paths
DEMULTI_DIR = "demultiplexing"
ALIGN_DIR = "alignment"
VECTOR_DIR = "vectoring"

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
