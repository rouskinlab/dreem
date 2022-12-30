from itertools import repeat
from typing import Any, Tuple, Dict, Callable, Iterable, Iterator


def args_for_apply(func: Callable,
                   iter_args: Iterable[Tuple[Any]],
                   iter_kwargs: Iterable[Dict[str, Any]]):
    """
    Convert a list where each item is a tuple of positional arguments
    and a list where each item is a dict of keyword arguments
    into one iterable that can be passed to starmap.

    Arguments
    ---------
    function: callable
        Function to be called once for each set of args and kwargs.
    
    iter_args: iterable[tuple[any]]
        Iterable of tuples of positional arguments. Each item in iter_args is a
        tuple of all positional arguments for the function.
    
    iter_kwargs: iterable[dict[str, any]]
        Iterable of dicts of keyword arguments. Each item in iter_kwargs is a
        dict of all keyword arguments for the function.

    Returns
    -------
    apply_args: iterator[callable, tuple[any], dict[str, any]]
        Arguments that can be passed to the function apply
    """
    apply_args = zip(repeat(func), iter_args, iter_kwargs)
    return apply_args


def apply(func: Callable, args: Tuple[Any], kwargs: Dict[str, Any]):
    """
    Call a function by unpacking positional and keyword arguments.

    Arguments
    ---------
    func: callable
        The function to call
    
    args: tuple[any]
        The positional arguments of the function
    
    kwargs: dict[str, any]
        The keyword arguments of the function
    
    Returns
    -------
    ret: any
        The return value of func called with *args and **kwargs
    """
    ret = func(*args, **kwargs)
    return ret


def starstarmap(starmapper: Callable[[Callable, Iterator], Iterable],
                func: Callable,
                iter_args: Iterable[Tuple[Any]],
                iter_kwargs: Iterable[Dict[str, Any]]):
    """
    Call a function multiple times, once for each group of positional and
    keyword arguments in iter_args and iter_kwargs, respectively.
    
    Arguments
    ---------
    starmapper: callable
        The engine performing the starmapping, such as
        itertools.starmap or multiprocessing.Pool.starmap
    
    func: callable
        The function called by starmapper with the arguments
    
    iter_args: iterable[tuple[any]]
        Iterable of tuples of positional arguments. Each item in iter_args is a
        tuple of all positional arguments for the function.
    
    iter_kwargs: iterable[dict[str, any]]
        Iterable of dicts of keyword arguments. Each item in iter_kwargs is a
        dict of all keyword arguments for the function.
    
    Returns
    -------
    results: iterable
        Iterable of return values of starmapper function
    """
    results = starmapper(apply, args_for_apply(func, iter_args, iter_kwargs))
    return results
