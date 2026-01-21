from dataclasses import dataclass
import logging
import os
import pathlib
from typing import Callable, Dict, Generic, List, Literal, Tuple, TypeVar
import numpy as np
import xarray as xr
import random

from shnitsel.core._api_info import internal
from shnitsel.core.typedefs import StateTypeSpecifier

# TODO: The `pathlib.Path` part of the Union gets mangled to `pathlib._local.Path`
# in the `write_shnitsel_file()` accessor when generating using
# Python 3.13; unfortunately, `pathlib._local.Path` doesn't appear to exist for
# earlier Python versions and causes an error on `import shnitsel.xarray`.
# Given that `isinstance(pathlib.Path(), os.PathLike)`, the truncated type alias
# might be adequate; if so, please remove this notice.
# PathOptionsType = str | os.PathLike | pathlib.Path
PathOptionsType = str | os.PathLike


@dataclass
class LoadingParameters:
    """Class to hold certain parameters required at loading time at various points in the reading/
    import pipeline
    """

    # A dict containing the information, which input observable has which unit. If not provided, the loader will guess the units either based on the default values of that simulator or the data in `path`
    input_units: Dict[str, str] | None = None
    # Flag to set how errors during loading are reported
    error_reporting: Literal['log', 'raise'] = 'log'

    # Optionally provide a dict of trajectory ids, mapping the (absolut) posix-paths of trajectories to ids or a function to map the path to an integer id
    trajectory_id: Dict[str, int] | Callable[[pathlib.Path], int] | None = None

    # Optionally provide a list of state types/multiplicities or a function to assign them to a dataset
    state_types: (
        StateTypeSpecifier
        | List[StateTypeSpecifier]
        | Callable[[xr.Dataset], xr.Dataset]
        | None
    ) = None

    # List of the names of states or a function to label them or None and let the trajectory loader make an educated guess
    state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None

    # Logger to use for log outputs
    logger: logging.Logger | None = None


@internal()
def make_uniform_path(
    path: PathOptionsType | None,
) -> pathlib.Path:
    """Unify the path options to alyways yield a pathlib.Path object

    Parameters
    ----------
    path : str | os.PathLike | pathlib.Path | None
        path input of arbitrary type

    Returns
    -------
    pathlib.Path|None
        The converted path

    Raises
    ------
    ValueError
        If path was `None` or could not be converted
    """
    if path is None:
        raise ValueError("Cannot canonize path `None`. Please provide a valid path.")

    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)
    return path


T = TypeVar("T", bound=int | float | complex | np.integer | np.floating)


class ConsistentValue(Generic[T]):
    """Class to keep track of a value that may only be assigned once and not overwritten afterwards.

    Can be used to check consistency of a value across multiple datasets.
    The value is written to and read from the property `v` of the object.

    Raises
    ------
    AttributeError
        Will be raised if the value is read before first assignment if the object has not been created with ``weak=true``.
    ValueError
        If an inconsistent value is assigned to this instance, i.e. two different values have been assigned.

    """

    def __init__(self, name="ConsistentValue", weak=False, ignore_none=False):
        self.name: str = name
        self.defined: bool = False
        self._weak: bool = weak
        self._val: T | None = None
        self._ignore_none: bool = ignore_none

    @property
    def v(self) -> T | None:
        if self.defined:
            return self._val
        elif self._weak:
            return None
        raise AttributeError(f"{self.name}.v accessed before assignment")

    @v.setter
    def v(self, new_val: T | None):
        if self._ignore_none and new_val is None:
            return

        if self.defined and new_val != self._val:
            raise ValueError(
                f"""inconsistent assignment to {self.name}:
current value: {type(self._val).__name__} = {repr(self._val)}
new value:  {type(new_val).__name__} = {repr(new_val)}
"""
            )

        self.defined = True
        self._val = new_val


def get_triangular(original_array: np.ndarray):
    """
    get_triangular - get the upper triangle of a (nstat1 x nstat2 x natoms x 3) matrix

    This function takes in a 4-dimensional numpy array (original_array) and returns a 3-dimensional numpy array (upper_tril)
    which is the upper triangle of the input matrix, obtained by excluding the diagonal elements.
    The number of steps (k) to move the diagonal above the leading diagonal is 1.
    The returned matrix has shape (len(cols), natoms, 3)

    Parameters
    ----------
    original_array : np.ndarray
        4D numpy array of shape (nstat1, nstat2, natoms, 3) representing the input matrix

    Returns
    -------
    upper_tril : np.ndarray
        3D numpy array of shape (len(cols), natoms, 3) representing the upper triangle of the input matrix
    """
    # Get the indices of the upper triangle
    nstat1, nstat2, natoms, xyz = original_array.shape

    if nstat1 != nstat2:
        raise ValueError("expected square input matrix")

    rows, cols = np.triu_indices(nstat2, k=1)
    upper_tril = np.zeros((len(cols), natoms, 3))

    for i in range(len(cols)):
        me = original_array[rows[i], cols[i]]
        upper_tril[i] = me

    return upper_tril


def dip_sep(dipoles: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Separates a complete matrix of dipoles into permanent
    and transitional dipoles, removing redundancy in the process.

    Parameters
    ----------
    dipoles : np.ndarray
        3D numpy array of shape (nstates, nstates, 3) where
        the first axis represents state before transition,
        the second axis represents state after transition and
        the third axis contains x, y and z coordinates.

    Returns
    -------
    dip_perm : np.ndarray
        2D numpy array of shape (nstates, 3)
    dip_trans : np.ndarray
        2D numpy array of shape (math.comb(nstates, 2), 3)
        in the order e.g. (for nstates = 4)
        0->1, 0->2, 0->3, 1->2, 1->3, 2->3
        where 0->1 is the transitional dipole between
        state 0 and state 1.
    """
    assert dipoles.ndim == 3
    nstates, check, three = dipoles.shape
    assert nstates == check
    assert three == 3
    dip_perm = np.diagonal(dipoles).T
    dip_trans = dipoles[np.triu_indices(nstates, k=1)]
    # logging.debug("permanent dipoles\n" + str(dip_perm))
    # logging.debug("transitional dipoles\n" + str(dip_trans))
    return dip_perm, dip_trans


def random_trajid_assigner(path: pathlib.Path) -> int:
    """Function to generate a random id for a path.

    Parameters
    ----------
    path : pathlib.Path
        Unused: the path we are generating for

    Returns
    -------
    int
        the chosen trajectory id
    """

    return random.randint(0, 2**31 - 1)
