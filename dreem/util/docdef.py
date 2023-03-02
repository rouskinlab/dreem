from functools import wraps
from inspect import getmembers, signature
import logging
from typing import Any, Callable

from click import Option

from ..util import cli


reserved_params = {"self", "cls"}


cli_options = dict(getmembers(cli, lambda member: isinstance(member, Option)))


cli_docs = {option.name: option.help for option in cli_options.values()}
api_docs = {
    "num_cpus": "Number of processors or threads to use",
    "fastq": "FASTQ file or pair of mated FASTQ files",
}
all_docs = {**api_docs, **cli_docs}
logging.debug(f"All automatic docstrings: {all_docs}")


api_defs = {
    "num_cpus": cli.NUM_CPUS,
}
# FIXME: Click options with no default value default to None. Options without a default value should not be added to cli_defs. But if an option is actually supposed to have a default value of None, then it will not be incorporated into cli_defs.
cli_defs = {option.name: option.default for option in cli_options.values()
            if option.default is not None}
all_defs = {**api_defs, **cli_defs}
logging.debug(f"All automatic definitions: {all_docs}")


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
        sig = signature(func)
        # Add information about every parameter to the docstring.
        param_lines = list()
        for name, param in sig.parameters.items():
            if name in reserved_params:
                # Ignore reserved parameters (if any).
                continue
            if doc := param_text.get(name):
                # Add */** for variable positional/keyword arguments.
                if param.kind == param.VAR_POSITIONAL:
                    name = f"*{name}"
                elif param.kind == param.VAR_KEYWORD:
                    name = f"**{name}"
                if param.annotation is not param.empty:
                    # Add the type annotation (if any) after the name
                    # of the parameter.
                    try:
                        name = f"{name}: {param.annotation.__name__}"
                    except AttributeError:
                        # Some types (e.g. UnionType) have no name.
                        name = f"{name}: {param.annotation}"
                if param.default is not param.empty:
                    # Add the default value (if any) in parentheses
                    # after the description of the parameter.
                    doc = f"{doc} (default: {repr(param.default)})"
                # Add the name and description of the parameter to the
                # docstring, on separate lines.
                param_lines.extend([f"{name}", f"    {doc}"])
            else:
                logging.warning("Missing documentation for parameter "
                                f"'{name}' of function '{func.__name__}'")
        # Assemble the new docstring.
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
        func.__doc__ = "\n".join(doc_lines)
        return func

    return decorator


def autodoc(ret_doc: str = "", **extras: str):
    """ Call ```paramdoc``` and automatically infer descriptions and
    type annotations about all parameters from the CLI and API. Extra
    annotations (if needed) may be given as keyword arguments. """
    return paramdoc(ret_doc, **all_docs, **extras)


def paramdef(*exclude, **defaults: Any):
    """ Give the keyword argments of a function default values. """

    def decorator(func: Callable):
        params = signature(func).parameters
        new_defaults = dict()
        for name, default in defaults.items():
            if (param := params.get(name)) is not None:
                if name in reserved_params:
                    logging.debug(f"Skipping reserved parameter '{name}' "
                                  f"in function '{func.__name__}'")
                    continue
                if name in exclude:
                    logging.debug(f"Skipping excluded parameter '{name}' "
                                  f"in function '{func.__name__}'")
                    continue
                if param.kind not in (param.POSITIONAL_OR_KEYWORD,
                                      param.KEYWORD_ONLY):
                    logging.debug(f"Skipping {param.kind} parameter '{name}' "
                                  f"in function '{func.__name__}'")
                    continue
                new_defaults[name] = default
                logging.debug(
                    f"Set default value of parameter '{name}' in function "
                    f"'{func.__name__}' to {defaults[name]}")

        @wraps(func)
        def new_func(*args, **kwargs):
            return func(*args, **{**new_defaults, **kwargs})

        return new_func

    return decorator


def autodef(*exclude, **extras):
    """ Call ```paramdef``` and automatically infer default values from
    the CLI and API. Extra defaults (if needed) may be given as keyword
    arguments. """
    return paramdef(*exclude, **all_defs, **extras)
