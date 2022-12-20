import os


# Global parameters
NUM_PROCESSES = cpus if (cpus := os.cpu_count()) else 1
BASE_COLORS = {"A": "#D3822A", "C": "#5AB4E5", "G": "#EAE958", "T": "#357766"}
PHRED_ENCODING = 33
BUFFER_LENGTH = 10000
