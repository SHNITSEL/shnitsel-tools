from typing import Any, overload

import numpy as np
import xarray as xr
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.geo.alignment import kabsch
from shnitsel.geo.analogs import extract_analogs
from shnitsel.analyze.generic import center

# def center_frames(da:xr.DataArray):
#     return da.groupby('frame').map(lambda x: center(x, dim='atom'))


@overload
def calc_rmsd(
    atXYZ: xr.DataArray | xr.Dataset | ShnitselDataset,
    reference: xr.DataArray,
    ignore_H: bool = True,
) -> xr.DataArray: ...


@overload
def calc_rmsd(
    atXYZ: TreeNode[Any, xr.DataArray | xr.Dataset | ShnitselDataset],
    reference: xr.DataArray,
    ignore_H: bool = True,
) -> TreeNode[Any, xr.DataArray]: ...


def calc_rmsd(
    atXYZ: xr.DataArray
    | xr.Dataset
    | ShnitselDataset
    | TreeNode[Any, xr.DataArray | xr.Dataset | ShnitselDataset],
    reference: xr.DataArray,
    ignore_H: bool = True,
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Function to calculate the Root Mean Squared Distance (RMSD) for a set of goemetries relative
    to a reference geometry.

    Parameters
    ----------
    atXYZ : xr.DataArray | xr.Dataset | ShnitselDataset | TreeNode[Any, xr.Dataset |xr.Dataset | ShnitselDataset]
        The source for geometries to calulate the RMSD for. If not provided as a data array, will be looked up in
        the variable `atXYZ` or `positions` of the dataset.
        If provided as a tree, the calc_rmsd will be mapped over the entries in the tree with the same reference.
    reference : xr.DataArray
        The geometry to be used as a reference. All dimensions of this reference geometry must be present in the
        geometries in `atXYZ`.
    ignore_H : bool, default=True
        To reduce the rmsd due to ambiguity in mapping multiple H atoms to their counterparts in the reference module,
        if set to `True`, the RMSD will be calculated only on non-H atoms.

    Returns
    -------
    xr.DataArray | TreeNode[Any, xr.DataArray]
        The RMSD result in the same (hierarchical) structure as the input.
    """

    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(calc_rmsd, reference=reference, ignore_H=ignore_H)

    geometry_data: xr.DataArray
    if isinstance(atXYZ, xr.DataArray):
        geometry_data = atXYZ
    else:
        dataset = wrap_dataset(atXYZ)
        geometry_data = dataset.positions

    geometry_data: xr.DataArray
    if not isinstance(reference, xr.DataArray):
        try:
            reference = reference.atXYZ
        except:
            logging.error(
                "Could not find geometry information in variable `atXYZ` in non-datarray parameter `reference`. Please provide reference geometry as a xr.DataArray."
            )
            raise

    assert {'atom', 'direction'}.issubset(geometry_data.dims), (
        "Geometries did not have `atom` and `direction` dimension."
    )
    assert {'atom', 'direction'}.issubset(reference.dims), (
        "Reference geometry did not have `atom` and `direction` dimension."
    )
    assert set(reference.dims).issubset(geometry_data.dims), (
        "Geometries for rmsd calculation did have all dimensions present on the reference geometry."
    )

    analog_geometries, analog_reference = extract_analogs([geometry_data, reference])

    aligned_analog_geometries = kabsch(
        analog_geometries, reference_or_indexers=analog_reference
    )

    aligned_distance = aligned_analog_geometries - analog_reference

    if ignore_H and 'atNums' in aligned_distance.coords:
        # Filter out H entries
        aligned_distance = aligned_distance.where(
            aligned_distance.atNums != 1, drop=True
        )

    rmsd: xr.DataArray = np.sqrt((aligned_distance**2).sum(['direction']).mean('atom'))

    return rmsd
