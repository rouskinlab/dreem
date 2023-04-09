from functools import wraps
from logging import getLogger
import os
from pathlib import Path
from shutil import rmtree
from typing import Callable


logger = getLogger(__name__)


# DREEM lock directory (for function lock_output)
LOCK_DIR = ".dreem-lock"


def lock_temp_dir(run: Callable):
    @wraps(run)
    def wrapper(*args, temp_dir: str | Path, save_temp: bool, **kwargs):
        lock_error = (f"The directory {temp_dir} is currently being used by "
                      f"another run of DREEM. If possible, use a different "
                      f"temporary directory that is not in use. If another "
                      f"run of DREEM crashed and left this directory locked, "
                      f"then you may maually delete {temp_dir}.")
        # Determine whether the temporary directory and the lock exist.
        lock = os.path.join(temp_dir, LOCK_DIR)
        try:
            os.mkdir(lock)
        except FileExistsError:
            # The lock already exists, which means another instance of
            # DREEM is using this temporary directory.
            raise SystemExit(lock_error)
        except FileNotFoundError:
            # The temporary directory does not exist yet, so create it
            # along with a lock.
            try:
                os.makedirs(lock, exist_ok=False)
            except FileExistsError:
                # If this error happens, it is due to a very unlikely
                # race condition wherein another run of DREEM raises a
                # FileNotFoundError from the step os.mkdir(lock), then
                # this run of DREEM does the same, then the first run
                # creates the directory with the step os.makedirs(lock),
                # and then this run tries to do the same thing but fails
                # because the directory was created moments before.
                raise SystemExit(lock_error)
            temp_dir_existed_before = False
            logger.debug(f"Created and locked temporary directory: {temp_dir}")
        else:
            # The temporary directory had existed, but the lock had not.
            temp_dir_existed_before = True
            logger.debug(f"Locked temporary directory: {temp_dir}")
        # The lock now exists, and any other run of DREEM that tires to
        # use the same lock will exit before it can use the temporary
        # directory or delete the lock. Thus, this run must delete the
        # lock when it exits, regardless of the circumstances.
        try:
            if temp_dir_existed_before and not save_temp:
                raise SystemExit(f"The directory {temp_dir} existed before "
                                 f"this run of DREEM was launched, and may "
                                 f"thus contain arbitrary files. Because DREEM "
                                 f"would delete this temporary directory after "
                                 f"finishing, causing unintentional loss of "
                                 f"any and all data in this directory, DREEM "
                                 f"is exiting as a precaution. Either specify "
                                 f"a temporary directory that does not exist "
                                 f"yet or use the option --save-temp (CLI) / "
                                 f"save_temp=True (API) to tell DREEM to not "
                                 f"delete the temporary directory on exiting.")
            # Run the wrapped function and capture its result.
            res = run(*args, temp_dir=temp_dir, save_temp=save_temp, **kwargs)
            if not save_temp:
                # If the run completes successfully and the temporary
                # directory should not be saved, then delete it.
                rmtree(temp_dir, ignore_errors=True)
                logger.debug(f"Deleted temporary directory: {temp_dir}")
            return res
        finally:
            # Always ensure that the temporary directory is unlocked
            # upon exiting.
            try:
                os.rmdir(lock)
                logger.debug(f"Unlocked temporary directory: {temp_dir}")
            except FileNotFoundError:
                pass
    # Return the decorator.
    return wrapper


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
