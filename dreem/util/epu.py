from typing import Any, Callable, Iterable, Sequence


def all_equal(values: Sequence[Any], same_types: bool = True):
    """ Return whether all values are equal. """
    if not values:
        # No values were given: return vacuous truth.
        return True
    val0 = values[0]
    if same_types:
        # Require that the data types match and the values are equal.
        typ0 = type(val0)
        return all(type(val) is typ0 and val == val0
                   for val in values[1:])
    else:
        # Require that the values are equal, but not the types.
        return all(val == val0 for val in values[1:])


def get_common_attrib(objects: Iterable[Any],
                      attrib: str, *,
                      call: bool = False,
                      key: Callable[[Any], Any] | None = None,
                      args: tuple[Any] | list[Any] = (),
                      kwargs: dict[str, Any] | None = None):
    """
    Return the consensus value of an attribute or method among one or
    more objects.

    Parameters
    ----------
    objects: Iterable[Any]
        One or more objects. Raise ValueError if empty.
    attrib: str
        Name of the attribute/method. Raise AttributeError if undefined.
    call: bool = False
        Whether the attribute is actually a method and should be called.
    key: Callable[[Any], Any] | None = None
        If given, check equality of the return values of this function,
        called with the attribute/method value as its single argument,
        rather than the attribute/method value itself.
    args: tuple[Any] | list[Any] = ()
        Positional arguments to pass to the method if call is True.
    kwargs: dict[str, Any] | None = None
        Keyword arguments to pass to the method if call is True.

    Returns
    -------
    Any
        If all EmClustering objects have the same value of the attribute
        or return the same value from the method, return that consensus
        value. Raise ValueError otherwise.
    """
    if kwargs is None:
        kwargs = dict()
    # Get the value of each object's attribute/method.
    values = [obj.__getattribute__(attrib) for obj in objects]
    if not values:
        raise ValueError("No objects were given")
    if call:
        # Call the attribute/method if call is True.
        values = [value(*args, **kwargs) for value in values]
    if key is None:
        # Compare the values.
        compare = values
    else:
        # Compare the return values of the key function.
        compare = list(map(key, values))
    if not all_equal(compare, same_types=True):
        raise ValueError(f"Got multiple values for '{attrib}': {values}")
    return values[0]
