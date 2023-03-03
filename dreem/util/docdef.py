from functools import wraps
from inspect import getmembers, Parameter, Signature
import logging
from shutil import get_terminal_size
from textwrap import dedent, wrap
from typing import Any, Callable

from click import Option

from ..util import cli

# Ignore special parameters with reserved names.
reserved_params = {"self", "cls"}

# Get every option defined for the command line interface.
cli_options = dict(getmembers(cli, lambda member: isinstance(member, Option)))

# Get the default value for every parameter.
api_defs = {
    "num_cpus": cli.NUM_CPUS,
}
cli_defs = {option.name: option.default for option in cli_options.values()
            if option.default is not None}
all_defs = {**cli_defs, **api_defs}

# Get the documentation for every parameter.
cli_docs = {option.name: option.help for option in cli_options.values()}
api_docs = {
    "kwargs": "Additional keyword-only arguments",
    "num_cpus": "Number of processors or threads to use",
    "fastq": "FASTQ file or pair of mated FASTQ files",
}
all_docs = {**cli_docs, **api_docs}


def get_param_default(param: Parameter,
                      defaults: dict[str, Any],
                      exclude_defs: tuple[str, ...]):
    """ Get a copy of the parameter, possibly with a new default. """
    if param.name in exclude_defs:
        logging.debug(f"Skipped excluded parameter '{param.name}'")
        return param
    if param.name in reserved_params:
        logging.debug(f"Skipped reserved parameter '{param.name}'")
        return param
    if param.kind != param.KEYWORD_ONLY:
        logging.debug(
            f"Skipped {param.kind.description} parameter '{param.name}'")
        return param
    try:
        default = defaults[param.name]
    except KeyError:
        logging.debug(f"Skipped parameter '{param.name}' with no default value")
        return param
    logging.debug(
        f"Set default value of parameter '{param.name}' to {repr(default)}")
    # Return copy of parameter with new default value.
    return param.replace(default=default)


def paramdef(defaults: dict[str, Any], exclude_defs: tuple[str, ...]):
    """ Give the keyword argments of a function default values. """
    if defaults is None:
        defaults = dict()

    def decorator(func: Callable):
        logging.debug("Setting default values of parameters for function "
                      f"'{func.__name__}'")
        # List all the parameters of the function, replacing the default
        # value of those parameters with defaults given in defaults.
        sig = Signature.from_callable(func)
        new_params = [get_param_default(param, defaults, exclude_defs)
                      for param in sig.parameters.values()]

        # Update the help text (does not affect actual default values).
        try:
            func.__signature__ = Signature(parameters=new_params)
        except ValueError as error:
            raise ValueError(f"Failed to set signature of {func.__name__} with "
                             f"parameters {', '.join(map(str, new_params))}. "
                             f"Raised error: {error}")

        # Update the actual default values (does not affect help text).
        new_defaults = {param.name: param.default for param in new_params}

        @wraps(func)
        def new_func(*args, **kwargs):
            return func(*args, **{**new_defaults, **kwargs})

        return new_func

    return decorator


def autodef(extra_defs: dict[str, Any] | None = None,
            exclude_defs: tuple[str, ...] = ()):
    """ Call ```paramdef``` and automatically infer default values from
    the CLI and API. Extra defaults (if needed) may be given as keyword
    arguments. """
    if extra_defs is None:
        extra_defs = dict()
    return paramdef({**all_defs, **extra_defs}, exclude_defs)


def get_param_lines(func: Callable, param_docs: dict[str, str]):
    sig = Signature.from_callable(func)
    # Add information about every parameter to the docstring.
    param_lines = list()
    for name, param in sig.parameters.items():
        if name in reserved_params:
            # Ignore reserved parameters (if any).
            continue
        if doc := param_docs.get(name):
            # Add the type annotation (if any) after the name of the
            # parameter.
            if param.annotation is param.empty:
                name_type = name
            else:
                try:
                    name_type = f"{name}: {param.annotation.__name__}"
                except AttributeError:
                    # Some types (e.g. UnionType) have no name.
                    name_type = f"{name}: {param.annotation}"
            # Add the default value (if any) in brackets after the
            # documentation of the parameter.
            if param.default is not param.empty:
                doc = f"{doc}  [default: {repr(param.default)}]"
            # Add the parameter's name, type, kind, and documentation to
            # the docstring.
            param_lines.extend([f"{name_type}  [{param.kind.description}]",
                                f"    {doc}"])
        else:
            logging.warning("Missing documentation for parameter "
                            f"'{name}' of function '{func.__name__}'")
    return param_lines


def get_doc_lines(func: Callable, param_lines: list[str], return_doc: str):
    sig = Signature.from_callable(func)
    doc_lines = list()
    if func.__doc__:
        # Use the existing docstring to start the new docstring.
        doc_lines.append(dedent(func.__doc__))
    else:
        logging.warning(f"Function '{func.__name__}' had no docstring")
    if param_lines:
        if doc_lines:
            doc_lines.append("")
        doc_lines.extend(["Parameters",
                          "----------"])
        doc_lines.extend(param_lines)
    if sig.return_annotation is not sig.empty:
        if return_doc:
            if doc_lines:
                doc_lines.append("")
            doc_lines.extend(["Return",
                              "------",
                              f"{sig.return_annotation}",
                              f"    {return_doc}"])
        else:
            logging.warning(
                f"Function '{func.__name__}' has return annotation "
                f"{sig.return_annotation} but was given no return text")
    elif return_doc:
        logging.warning(f"Function '{func.__name__}' was given return text "
                        f"'{return_doc}' but has no return annotation")
    return doc_lines


def paramdoc(param_docs: dict[str, str], return_doc: str, width: int):
    """
    Give a function a new docstring where each parameter gets annotated
    with the text given in the keyword arguments (```param_docs```).

    Parameters
    ----------
    param_docs: dict[str, str]
        Description of each parameter, keyed by name
    return_doc: str
        Description of the return value; will occur at end of docstring
    width: int
        Length to which to wrap the text, or 0 for no wrapping

    Return
    ------
    callable
        Function with new docstring
    """

    def decorator(func: Callable):
        param_lines = get_param_lines(func, param_docs)
        doc_lines = get_doc_lines(func, param_lines, return_doc)
        doc_str = "\n".join(doc_lines)
        func.__doc__ = wrap(doc_str, width=width) if width >= 1 else doc_str
        return func

    return decorator


def autodoc(extra_docs: dict[str, str] | None = None, return_doc: str = ""):
    """ Call ```paramdoc``` and automatically infer descriptions and
    type annotations about all parameters from the CLI and API.
    Documentation of any extra parameters may also be given. """
    if extra_docs is None:
        extra_docs = dict()
    return paramdoc({**all_docs, **extra_docs},
                    return_doc,
                    get_terminal_size().columns)


def auto(*,
         extra_defs: dict[str, Any] | None = None,
         exclude_defs: tuple[str, ...] = (),
         extra_docs: dict[str, str] | None = None,
         return_doc: str = ""):
    """ Combine ```autodef``` and ```autodoc```, in that order. """
    def decorator(func: Callable):
        return autodoc(extra_docs, return_doc)(
            autodef(extra_defs, exclude_defs)(func)
        )

    return decorator
