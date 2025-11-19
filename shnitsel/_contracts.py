from typing import Literal, Callable
from collections import namedtuple

_fields = [
    'to_be',
    'dims',
    'coords',
    'data_vars',
    'groupable',
    'coords_or_vars',
    'name',
    'not_dims',
]
Needs = namedtuple('Needs', _fields, defaults=(None,) * len(_fields))

# TODO: FIXME: This decorator breaks type hints and autocomplete on my machine. I still see the documentation, but it no longer suggests parameter names.


def needs(
    to_be: Literal['da', 'ds', None] = None,
    dims: set[str] | None = None,
    coords: set[str] | None = None,
    data_vars: set[str] | None = None,
    groupable: set[str] | None = None,
    coords_or_vars: set[str] | None = None,
    name: str | None = None,
    not_dims: set[str] | None = None,
) -> Callable[[Callable], Callable]:
    def decorator(func: Callable) -> Callable:
        if to_be == 'da' and data_vars is not None:
            raise ValueError()
        func._needs = Needs(
            to_be, dims, coords, data_vars, groupable, coords_or_vars, name, not_dims
        )
        return func

    return decorator
