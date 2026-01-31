from typing import Concatenate, Literal, Callable, ParamSpec, TypeVar
from collections import namedtuple
import xarray as xr

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

T = TypeVar('T')
P = ParamSpec('P')

DatasetOrArray = TypeVar('DatasetOrArray', bound=xr.Dataset | xr.DataArray)


def needs(
    to_be: Literal['da', 'ds', None] = None,
    dims: set[str] | None = None,
    coords: set[str] | None = None,
    data_vars: set[str] | None = None,
    groupable: set[str] | None = None,
    coords_or_vars: set[str] | None = None,
    name: str | None = None,
    not_dims: set[str] | None = None,
) -> Callable[
    [Callable[Concatenate[DatasetOrArray, P], T]],
    Callable[Concatenate[DatasetOrArray, P], T],
]:
    """
    Decorator to add information about required data in an initial argument of a function that can be applied to a Dataset 
    or a DataArray.

    Parameters
    ----------
    to_be : Literal['da', 'ds', None] , optional
            Whether the initial argument needs to be a DataArray or a Dataset. Defaults to None.
    dims : set[str] | None, optional
            A list of dimensions that need to be present. Defaults to None.
    coords : set[str] | None, optional
            A list of coordinates that need to be present. Defaults to None.
    data_vars : set[str] | None, optional
            A list of data variables that need to be present. Defaults to None.
    groupable : set[str] | None, optional
            A list of coordinates by which the data needs to be groupable. Defaults to None.
    coords_or_vars : set[str] | None, optional
            A list of entries that can either be present as a coordinate or a variable. Defaults to None.
    name : str | None, optional
            Whether the data needs to have an attribute `name` set. Defaults to None.
    not_dims : set[str] | None, optional
            A list of dimensions that must not be present on the data. Defaults to None.

    Returns
    -------
    Callable[[Callable[P,T]], Callable[P,T]]
        A descriptor that does not change the method signature but sets a `_needs` property of the function object for later use.
    """

    def decorator(
        func: Callable[Concatenate[DatasetOrArray, P], T],
    ) -> Callable[Concatenate[DatasetOrArray, P], T]:
        if to_be == 'da' and data_vars is not None:
            raise ValueError(
                "If a DataArray is required, you cannot specify `data_vars` as requirements as there is only one Variable in each DataArray."
            )

        func._needs = Needs(
            to_be, dims, coords, data_vars, groupable, coords_or_vars, name, not_dims
        )
        return func

    return decorator
