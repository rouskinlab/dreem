from functools import wraps
from inspect import getmembers, Parameter, Signature
import logging
from typing import Any, Callable

from click import Option

from ..util import cli

reserved_params = {"self", "cls"}

cli_options = dict(getmembers(cli, lambda member: isinstance(member, Option)))

cli_docs = {option.name: option.help for option in cli_options.values()}
api_docs = {
    "kwargs": "Additional keyword-only arguments",
    "num_cpus": "Number of processors or threads to use",
    "fastq": "FASTQ file or pair of mated FASTQ files",
}
all_docs = {**api_docs, **cli_docs}
logging.debug(f"All automatic docstrings: {all_docs}")

api_defs = {
    "num_cpus": cli.NUM_CPUS,
}
cli_defs = {option.name: option.default for option in cli_options.values()
            if option.default is not None}
all_defs = {**api_defs, **cli_defs}
logging.debug(f"All automatic definitions: {all_docs}")


def get_param_lines(func: Callable, **param_text: str):
    sig = Signature.from_callable(func)
    # Add information about every parameter to the docstring.
    param_lines = list()
    for name, param in sig.parameters.items():
        if name in reserved_params:
            # Ignore reserved parameters (if any).
            continue
        if doc := param_text.get(name):
            if param.annotation is param.empty:
                title = f"{name}:"
            else:
                # Add the type annotation (if any) after the name
                # of the parameter.
                try:
                    title = f"{name}: {param.annotation.__name__}"
                except AttributeError:
                    # Some types (e.g. UnionType) have no name.
                    title = f"{name}: {param.annotation}"
            if param.default is not param.empty:
                # Add the default value (if any) in parentheses
                # after the description of the parameter.
                doc = f"{doc} (default: {repr(param.default)})"
            # Add the name and description of the parameter to the
            # docstring, on separate lines.
            param_lines.extend([f"{title} ({param.kind.description})",
                                f"    {doc}"])
        else:
            logging.warning("Missing documentation for parameter "
                            f"'{name}' of function '{func.__name__}'")
    return param_lines


def get_doc_lines(func: Callable, param_lines: list[str], ret_desc: str):
    sig = Signature.from_callable(func)
    doc_lines = list()
    if func.__doc__:
        # Use the existing docstring to start the new docstring.
        doc_lines.append(func.__doc__.strip())
    else:
        logging.warning(f"Function '{func.__name__}' had no docstring")
    if param_lines:
        if doc_lines:
            doc_lines.append("")
        doc_lines.extend(["Parameters",
                          "----------"])
        doc_lines.extend(param_lines)
    if sig.return_annotation is not sig.empty:
        if ret_desc:
            if doc_lines:
                doc_lines.append("")
            doc_lines.extend(["Return",
                              "------",
                              f"{sig.return_annotation}",
                              f"    {ret_desc}"])
        else:
            logging.warning(
                f"Function '{func.__name__}' has return annotation "
                f"{sig.return_annotation} but was given no return text")
    elif ret_desc:
        logging.warning(f"Function '{func.__name__}' was given return text "
                        f"'{ret_desc}' but has no return annotation")
    return doc_lines


def paramdoc(ret_desc: str = "", **param_text: str):
    """
    Give a function a new docstring where each parameter gets annotated
    with the text given in the keyword arguments (```param_text```).

    Parameters
    ----------
    ret_desc: str
        Description of the return value; will occur at end of docstring
    **param_text: str
        Description of each parameter, keyed by name

    Return
    ------
    callable
        Function with new docstring
    """

    def decorator(func: Callable):
        param_lines = get_param_lines(func, **param_text)
        doc_lines = get_doc_lines(func, param_lines, ret_desc)
        func.__doc__ = "\n".join(doc_lines)
        return func

    return decorator


def autodoc(ret_doc: str = "", **extras: str):
    """ Call ```paramdoc``` and automatically infer descriptions and
    type annotations about all parameters from the CLI and API. Extra
    annotations (if needed) may be given as keyword arguments. """
    return paramdoc(ret_doc, **all_docs, **extras)


def get_param_default(param: Parameter,
                      defaults: dict[str, Any],
                      exclude: tuple[str, ...]):
    """ Get a copy of the parameter, possibly with a new default. """
    if param.name in exclude:
        logging.debug(f"Skipped excluded parameter '{param.name}'")
        return param.replace()
    if param.name in reserved_params:
        logging.debug(f"Skipped reserved parameter '{param.name}'")
        return param.replace()
    if param.kind != param.KEYWORD_ONLY:
        logging.debug(
            f"Skipped {param.kind.description} parameter '{param.name}'")
        return param.replace()
    try:
        default = defaults[param.name]
    except KeyError:
        logging.debug(f"Skipped parameter '{param.name}' with no default value")
        return param.replace()
    logging.debug(
        f"Set default value of parameter '{param.name}' to {repr(default)}")
    # Return copy of parameter with new default value.
    return param.replace(default=default)


def paramdef(defaults: dict[str, Any] | None = None,
             exclude: tuple[str, ...] = ()):
    """ Give the keyword argments of a function default values. """
    if defaults is None:
        defaults = dict()

    def decorator(func: Callable):
        logging.debug("Setting default values of parameters for function "
                      f"'{func.__name__}'")
        # List all the parameters of the function, replacing the default
        # value of those parameters with defaults given in defaults.
        sig = Signature.from_callable(func)
        new_params = [get_param_default(param, defaults, exclude)
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
            return func(*args, {**new_defaults, **kwargs})

        return new_func

    return decorator


def autodef(extras: dict[str, Any] | None = None,
            exclude: tuple[str, ...] = ()):
    """ Call ```paramdef``` and automatically infer default values from
    the CLI and API. Extra defaults (if needed) may be given as keyword
    arguments. """
    if extras is None:
        extras = dict()
    return paramdef({**all_defs, **extras}, exclude)
