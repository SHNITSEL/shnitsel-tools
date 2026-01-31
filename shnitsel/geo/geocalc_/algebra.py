import xarray as xr

from shnitsel.core.typedefs import DataArrayOrVar
from shnitsel.analyze.generic import norm
import numpy as np


def normalize(a: DataArrayOrVar, dim='direction') -> DataArrayOrVar:
    """Return vectors normalized along the `dim` dimension.

    Parameters
    ----------
    a : xr.DataArray | xr.Variable
        The data array to normalize along the `dim` dimension.

    Returns
    -------
    xr.DataArray | xr.Variable
        Resulting DataArray after calculation of the noramlized vectors in the dimension `dim`.
    """
    return a / norm(a, dim=dim)


def dnorm(a: DataArrayOrVar) -> DataArrayOrVar:
    """Calculate the norm along the `direction` dimension. All other dimensions are maintaned

    Parameters
    ----------
    a : xr.DataArray | xr.Variable
        The data array to perform the norming on

    Returns
    -------
    xr.DataArray | xr.Variable 
        Resulting dataarray after calculation of the norm in the dimension `direction.
    """
    return norm(a, dim='direction')


def dcross(a: DataArrayOrVar, b: DataArrayOrVar) -> DataArrayOrVar:
    """Generalized cross vector product in the dimension of `direction`.

    Parameters
    ----------
    a : xr.DataArray | xr.Variable
        The first array to use for the binary operation
    b : xr.DataArray | xr.Variable
        The second array to use for the binary operation

    Returns
    -------
    xr.DataArray | xr.Variable
        The resulting array of the cross-product
    """
    return xr.cross(a, b, dim='direction')  # type: ignore # Type should be same as input


def ddot(a: DataArrayOrVar, b: DataArrayOrVar) -> DataArrayOrVar:
    """Dot product in the dimension of `direction`.

    Parameters
    ----------
    a : xr.DataArray | xr.Variable 
        The first array to use for the binary operation
    b : xr.DataArray | xr.Variable 
        The second array to use for the binary operation

    Returns
    -------
    xr.DataArray | xr.Variable
        The resulting array of the dot-product still retaining all other dimensions except `direction`.
    """
    return xr.dot(a, b, dim='direction')


def angle_(a: DataArrayOrVar, b: DataArrayOrVar) -> DataArrayOrVar:
    """Helper function to calculate the angle between the entries in a and b based on their coordinates in the `direction` dimension.

    Parameters
    ----------
    a : xr.DataArray | xr.Variable
        The first array to use for the binary operation
    b : xr.DataArray | xr.Variable
        The second array to use for the binary operation

    Returns
    -------
    xr.DataArray | xr.Variable
        The resulting array of the angle calculation still retaining all other dimensions except `direction`.
    """

    cos_theta = ddot(a, b) / (dnorm(a) * dnorm(b))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.arccos(cos_theta)


def angle_cos_sin_(
    a: DataArrayOrVar, b: DataArrayOrVar
) -> tuple[DataArrayOrVar, DataArrayOrVar]:
    """
    Returns the cosine and sine of the angle between two vectors

    Parameters
    ----------
    a : DataArrayOrVar
        First vector
    b : DataArrayOrVar
        Second vector

    Returns
    -------
    tuple[DataArrayOrVar, DataArrayOrVar]
        First the vector of cosines, second the vector of sines.
    """
    prod = dnorm(a) * dnorm(b)
    return (
        ddot(a, b) / prod,
        dnorm(dcross(a, b)) / prod,
    )


def normal(
    a: DataArrayOrVar,
    b: DataArrayOrVar,
    c: DataArrayOrVar,
) -> DataArrayOrVar:
    """Calculate normal vectors on the planes through corresponding points in a, b and c.

    The normal vector will be calculated based on the position in `direction` dimension.

    Parameters
    ----------
    a : xr.DataArray | xr.Variable
        The first array to use for the ternary operation
    b : xr.DataArray | xr.Variable
        The second array to use for the ternary operation
    c : xr.DataArray | xr.Variable
        The third array to use for the ternary operation

    Returns
    -------
    xr.DataArray | xr.Variable
        An array with all dimensions equal to those of a, b, and c but holding normal vectors along the `direction` dimension.
    """
    return dcross(a - b, c - b)
