from typing import Iterable, Literal, Sequence
from shnitsel._contracts import needs
import xarray as xr
import numpy as np

from shnitsel.core._api_info import API
from shnitsel.core.typedefs import AtXYZ


@API()
@needs(dims={'atom'})
def get_centered_geometry(atXYZ: AtXYZ, by_mass: Literal[False] = False) -> AtXYZ:
    """Helper function to set the center of the geometry (i.e. mean along the `atom` axis) to zero.

    Parameters
    ----------
    atXYZ : AtXYZ
        Array of positional data
    by_mass : Literal[False], optional
        Flag whether the centering/average should be center of mass or just plain average of positions. Defaults to False.

    Raises
    ------
    NotImplementedError
        Centering the COM instead of the mean is currently not implemented.

    Returns
    -------
    AtXYZ
        Resulting positions after centering.
    """
    if by_mass:
        raise NotImplementedError
    return atXYZ - atXYZ.mean('atom')


def rotational_procrustes_np(
    A: np.ndarray, B: np.ndarray, weight: Sequence[float] | None = None
) -> np.ndarray:
    """Rotationally align the geometrie(s) in A to the single geometry in B.

    This helper function is specifically tailored to work directly on numpy arrays
    in contrast to its namesake `rotational_procrustes()`, which accepts
    `xarray.DataArray` parameters.

    Parameters
    ----------
    A : np.ndarray
        The geometries to process with shape
        ``(n_geometries, n_points, n_coordinates)``
    B : np.ndarray
        The reference geometry with shape
        ``(n_points, n_coordinates)``
    weight : Sequence[float], optional
        How much importance should be given to the alignment of
        each point, by default equal importance

    Returns
    -------
    np.ndarray
        An array with the same shape as A

    Notes
    -----
    We are solving the following minimization problem for each geometry `A(i)`:
    .. math::
        \\min_R \\Vert A(i)R - B\\Vert^2_F

    For this, we calculte the weighted cross-covariance matrix:
    We are solving the minimization problem
    .. math::
        C = A^T B

    And then we find the SVD
    .. math::
        C = U \\Sigma V^T

    Which provides us with the optimum rotation:
    .. math::
        R = V U^T

    If the resulting R has negative determinant, the rotation is instead
    a mirroring operation and we invert the sign of the last column of vt
    to restore rotational properties.
    """
    from scipy.linalg import svd

    if weight is not None:
        A = np.diag(weight) @ A

    # np.matrix_transpose always swaps last two axes, whereas
    # NDArray.T reverses the order of all axes.
    t = np.matrix_transpose
    # The following uses a double transpose in imitation of
    # scipy's orthogonal_procrustes, where this is said to
    # save memory. t(t(B) @ A) == t(A) @ B.
    u, _, vt = svd(t(t(B) @ A))
    # Flip the sign of the last row of each stacked vt matrix
    # depending on the sign of the corresponding determinant.
    # This is an alternative implementation of the algorithm
    # used in qcdev's procrustes.rotation.
    vt[..., -1, :] *= np.sign(np.linalg.det(u @ vt))[:, None]
    R = u @ vt
    return A @ R


def rotational_procrustes(
    A: xr.DataArray,
    B: xr.DataArray,
    dim0: str = 'atom',
    dim1: str = 'direction',
    weight: Sequence[float] | None = None,
) -> xr.DataArray:
    """Rotationally align the geometry or geometries in A to the single geometry in B.

    Parameters
    ----------
    A : xr.DataArray
        The (optionally multiple) geometries to process
    B : xr.DataArray
        The reference geometry
    dim0 : str, optional
        The name of the dimension over points to be rotated;
        must be present in ``A`` and ``B`; by default 'atom'
    dim1 : str, optional
        The name of the dimension over the coordinates of the aforementioned
        points; must be present in ``A`` and ``B`; by default 'direction'
    weight : Sequence[float], optional
        How much importance should be given to the alignment of
        each point (atom), by default equal importance

    Returns
    -------
    xr.DataArray
        An xr.DataArray with the same shape as ``A`` but with entries aligned to the overall geometry of ``B``
    """
    return xr.apply_ufunc(
        rotational_procrustes_np,
        A,
        B,
        input_core_dims=[[dim0, dim1], [dim0, dim1]],
        output_core_dims=[[dim0, dim1]],
        kwargs={'weight': weight},
    )


@needs(dims={'atom', 'direction'})
def kabsch(
    atXYZ: xr.DataArray,
    reference_or_indexers: xr.DataArray | dict | None = None,
    **indexers_kwargs,
) -> xr.DataArray:
    """Rotationally align the molecular geometries in ``atXYZ`` to a single molecular geometry.

    If no `reference_or_indexers` argument (or the `indexers_kwargs` option) is passed,
    this function will try to use the first frame or first timestep in `atXYZ` as a reference.

    Parameters
    ----------
    atXYZ : xr.DataArray
        The geometries to process (with dims 'atom', 'direction')

    reference_or_indexers : xr.DataArray | dict, optional
        Either a reference geometry (with dims 'atom', 'direction')
        or an indexer dictionary which will be passed to ``atXYZ.sel()``
        to indetify a single geometry in the `atXYZ` parameter to
        use as a reference point.

    **indexer_kwargs
        The keyword-argument form of the indexer to be passed to ``atXYZ.sel()``

    Returns
    -------
    xr.DataArray
        The aligned geometries

    Raises
    ------
    ValueError
        If nothing is done to indicate a reference geometry, i.e.
        neither reference_or_indexers nor indexer_kwargs are passed
    """

    if isinstance(reference_or_indexers, xr.DataArray):
        reference = reference_or_indexers
    elif isinstance(reference_or_indexers, dict):
        reference = atXYZ.sel(reference_or_indexers)
    elif len(indexers_kwargs) != 0:
        reference = atXYZ.sel(indexers_kwargs)
    elif 'frame' in atXYZ.dims:
        reference = atXYZ.isel(frame=0)
    elif 'time' in atXYZ.dims:
        reference = atXYZ.isel(time=0)
    else:
        raise ValueError(
            "Could not retrieve reference geometry. Please specify a reference geometry."
        )

    # TODO: is it ever necessary to center the molecule?
    # If so, should this always be done using the physical center-of-mass,
    # or is it ever appropriate to use the unweighted mean of points?

    # TODO: Kevin: As far as I understand the theory, you first center all geometries
    # Then you remove all scales via the mean RMSD scale and then you apply rotational procustes.
    # Based on the SVD, the S matrix should be able to figure out the scale, but we do not
    # Use that for rotation. May be that the centering is the key part that is missing.

    return rotational_procrustes(atXYZ, reference)
