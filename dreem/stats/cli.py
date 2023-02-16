from dreem.util.cli import *
from dreem.stats.main import run


@opt_out_dir
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == "__main__":
    cli()
