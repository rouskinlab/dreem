from functools import wraps
from logging import getLogger
from os import makedirs, path, rmdir
from typing import Callable


logger = getLogger(__name__)

# DREEM lock directory (for function lock_output)
LOCK_DIR = ".dreem-lock"


def lock_output(run: Callable):
    """ Prevent multiple instances of DREEM from using the same output
    and temporary directories simultaneously, to avoid data races. """
    @wraps(run)
    def with_locked_output(*args, out_dir, temp_dir, **kwargs):
        # Create hidden lock directories inside output and temp that
        # exist only for the duration of the run function.
        out_lock = path.join(out_dir, LOCK_DIR)
        temp_lock = path.join(temp_dir, LOCK_DIR)
        all_locks = [out_lock, temp_lock]
        my_locks = list()
        try:
            for lock in all_locks:
                try:
                    # Creating the locks will fail if another instance of
                    # DREEM (with its own locks) is currently running.
                    makedirs(lock, exist_ok=False)
                except FileExistsError:
                    # Quit because another instance of DREEM is running.
                    raise SystemExit(f"An instance of DREEM with directory "
                                     f"{path.dirname(lock)} is running, "
                                     f"which can cause data races.")
                else:
                    my_locks.append(lock)
                    logger.debug(f"Created directory lock: {lock}")
            try:
                # Call the DREEM run function and capture its return value.
                result = run(*args, out_dir=out_dir, temp_dir=temp_dir, **kwargs)
            except Exception:
                raise
            else:
                return result
        finally:
            # Always delete the locks after the run finishes, whether
            # normally or with an error.
            while my_locks:
                lock = my_locks.pop()
                try:
                    rmdir(lock)
                    logger.debug(f"Removed directory lock: {lock}")
                except FileNotFoundError:
                    logger.error(f"Failed to delete lock {lock}; please delete "
                                 f"this directory yourself, if it exists")
    # Return the decorator.
    return with_locked_output


def get_num_parallel(n_tasks: int,
                     max_procs: int,
                     parallel: bool,
                     hybrid: bool = False) -> tuple[int, int]:
    """ Determine how to parallelize the tasks.

    Parameters
    ----------
    n_tasks: int (≥ 1)
        Number of tasks to parallelize
    max_procs: int (≥ 1)
        Maximum number of processes to run at one time
    parallel: bool
        Whether to permit multiple tasks to be run in parallel
    hybrid: bool (default: False)
        Whether to allow both multiple tasks to run in parallel and,
        at the same, each task to run multiple processes in parallel

    Returns
    -------
    int (≥ 1)
        Number of tasks to run in parallel
    int (≥ 1)
        Number of processes to run for each task
    """
    if n_tasks >= 1 and max_procs >= 1:
        # This function only works if there is at least one task to
        # parallelize, at least one process is allowed, and parallel
        # is a valid option.
        if parallel:
            # Multiple tasks may be run in parallel. The number of tasks
            # run in parallel cannot exceed 1) the total number of tasks
            # and 2) the user-specified maximum number of processes.
            n_tasks_parallel = min(n_tasks, max_procs)
        else:
            # Otherwise, only one task at a time can be run.
            n_tasks_parallel = 1
        if n_tasks_parallel == 1 or hybrid:
            # Each individual task can be run by multiple processes in
            # parallel, as long as either 1) multiple tasks are not run
            # simultaneously in parallel (i.e. n_tasks_parallel == 1)
            # or 2) the calling function sets hybrid=True, which lets
            # multiple tasks run in parallel and each run with multiple
            # processes. Only the alignment module can simultaneously
            # run multiple tasks and multiple processes for each task
            # because its two most computation-heavy processes (cutadapt
            # and bowtie2) come with their own parallelization abilities
            # that can work independently of Python's multiprocessing
            # module. However, the other modules (e.g. vectoring) are
            # parallelized using the multiprocessing module, which does
            # not support "nesting" parallelization in multiple layers.
            # Because n_tasks_parallel is either 1 or the smallest of
            # n_tasks and n_procs (both of which are ≥ 1), it must be
            # that 1 ≤ n_tasks_parallel ≤ n_procs, and therefore that
            # 1 ≤ n_procs / n_tasks_parallel ≤ n_procs, so the
            # integer quotient must be a valid number of processes.
            n_procs_per_task = max_procs // n_tasks_parallel
        else:
            # Otherwise, only one process can work on each task.
            n_procs_per_task = 1
    else:
        logger.warning("Defaulting to 1 process due to invalid number of "
                       f"tasks ({n_tasks}) and/or processes ({max_procs}).")
        n_tasks_parallel = 1
        n_procs_per_task = 1
    return n_tasks_parallel, n_procs_per_task
