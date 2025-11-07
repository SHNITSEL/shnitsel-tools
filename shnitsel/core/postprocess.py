import math
import itertools
from logging import warning

from typing import Collection, Hashable, TypeAlias

import numpy as np
import xarray as xr

from .._contracts import needs
from . import xrhelpers

##############################################
# For backward compatibility while refactoring
from .convenience import (  # Consider these for deprecation
    pairwise_dists_pca as pairwise_dists_pca,
    pca_and_hops as pca_and_hops,
    relativize as relativize,
    hop_indices as hop_indices,
    ts_to_time as ts_to_time,
    setup_frames as setup_frames,
    validate as validate,
)
from .geom import (
    distance as distance,
    angle as angle,
    dihedral as dihedral,
    normal as normal,
)
from .populations import classical as calc_pops  # noqa: F401
from .ml import pca as pca
from .spectra import (
    assign_fosc,
    get_fosc as get_fosc,
    broaden_gauss as broaden_gauss,
    ds_broaden_gauss as ds_broaden_gauss,
)
from .stats import (
    calc_ci as calc_ci,
    ci_agg_last_dim as ci_agg_last_dim,
    xr_calc_ci as xr_calc_ci,
    time_grouped_ci as time_grouped_ci,
)

# Question marks for those that don't actually accept a mol object
from ..rd import (
    set_atom_props as set_atom_props,
    to_mol as to_mol,  # ?
    mol_to_numbered_smiles as mol_to_numbered_smiles,
    numbered_smiles_to_mol as numbered_smiles_to_mol,  # ?
    smiles_map as smiles_map,  # ?
    default_mol as default_mol,  # ?
)

# End of backward-compatability imports
#######################################

Astates: TypeAlias = xr.DataArray
AtXYZ: TypeAlias = xr.DataArray
DimName: TypeAlias = Hashable
Frames: TypeAlias = xr.Dataset
PerState: TypeAlias = xr.Dataset
InterState: TypeAlias = xr.Dataset

_var_delta_t_msg = "`delta_t` varies between the trajectories. Please separate the trajectories into groups"


def norm(
    da: xr.DataArray, dim: DimName = 'direction', keep_attrs: bool | str | None = None
) -> xr.DataArray:
    """Calculate the 2-norm of a DataArray, reducing the dimension with name *dim*

    Parameters
    ----------
    da
        Array to calculate the norm of
    dim, optional
        Dimension to calculate norm along (and therby reduce), by default 'direction'
    keep_attrs, optional
        How to deal with attributes; passed to xr.apply_ufunc, by default None

    Returns
    -------
        A DataArray with dimension *dim* reduced
    """
    res: xr.DataArray = xr.apply_ufunc(
        np.linalg.norm,
        da,
        input_core_dims=[[dim]],
        on_missing_core_dim='copy',
        kwargs={"axis": -1},
        keep_attrs=keep_attrs,
    )
    return res


def subtract_combinations(
    da: xr.DataArray, dim: DimName, labels: bool = False
) -> xr.DataArray:
    """Calculate all possible pairwise differences over a given dimension

    Parameters
    ----------
    da
        Input DataArray; must contain dimension `dim`
    dim
        Dimension (of size $n$) to take pairwise differences over
    labels, optional
        If True, label the pairwise differences based on the index of `dim`, by default False

    Returns
    -------
        A DataArray with the dimension `dim` replaced by a dimension '`dim`comb' of size $n(n-1)/2$
    """

    def midx(da, dim):
        return xrhelpers.midx_combs(da.get_index(dim))[f'{dim}comb']

    if dim not in da.dims:
        raise ValueError(f"'{dim}' is not a dimension of the DataArray")
    
    combination_dimension_name = f"{dim}comb"
    if combination_dimension_name in da:
        #TODO: FIXME: Appropriately remove the 'combination_dimension_name' indices and coordinates from the Dataset
        raise ValueError(f"'{combination_dimension_name}' is already an index, a variable or a coordinate of the DataArray")

    n = da.sizes[dim]

    mat = np.zeros((math.comb(n, 2), n))
    combs = itertools.combinations(range(n), 2)

    # After matrix multiplication, index r of output vector has value c2 - c1
    for r, (c1, c2) in enumerate(combs):
        mat[r, c1] = -1
        mat[r, c2] = 1

    if labels:
        xrmat = xr.DataArray(
            data=mat, coords={f'{dim}comb': midx(da, dim), dim: da.get_index(dim)}
        )
    else:
        xrmat = xr.DataArray(
            data=mat,
            dims=[f'{dim}comb', dim],
        )

    newdims = list(da.dims)
    newdims[newdims.index(dim)] = f'{dim}comb'

    res = (xrmat @ da).transpose(*newdims)
    res.attrs = da.attrs
    res.attrs['deltaed'] = set(res.attrs.get('deltaed', [])).union({dim})
    return res


@needs(dims={'frame'})
def sudi(da: xr.DataArray) -> xr.DataArray:
    """Take successive differences along the 'frame' dimension

    Parameters
    ----------
    da
        An ``xarray.DataArray`` with a 'frame' dimension corresponding
        to a ``pandas.MultiIndex`` of which the innermost level is 'time'.

    Returns
    -------
        An ``xarray.DataArray`` with the same shape, dimension names etc.,
        but with the data of the (i)th frame replaced by the difference between
        the original (i+1)th and (i)th frames, with zeros filling in for both the
        initial frame and any frame for which time = 0, to avoid taking differences
        between the last and first frames of successive trajectories.
    """
    res = xr.apply_ufunc(
        lambda arr: np.diff(arr, prepend=np.array(arr[..., [0]], ndmin=arr.ndim)),
        da,
        input_core_dims=[['frame']],
        output_core_dims=[['frame']],
    )
    res[{'frame': res['time'] == 0}] = 0
    return res


##############################################
# Functions generally applicable to timeplots:

def keep_norming(
    da: xr.DataArray, exclude: Collection[DimName] | None = None
) -> xr.DataArray:
    if exclude is None:
        exclude = {'state', 'statecomb', 'frame'}
    for dim in set(da.dims).difference(exclude):
        da = norm(da, dim, keep_attrs=True)
        da.attrs['norm_order'] = 2
    return da


@needs(dims={'state'})
def get_per_state(frames: Frames) -> PerState:
    props_per = {'energy', 'forces', 'dip_perm'}.intersection(frames.keys())
    per_state = frames[props_per].map(keep_norming, keep_attrs=False)
    per_state['forces'] = per_state['forces'].where(per_state['forces'] != 0)

    per_state['energy'].attrs['long_name'] = r'$E$'
    per_state['forces'].attrs['long_name'] = r'$\mathbf{F}$'
    if 'dip_perm' in per_state:
        per_state['dip_perm'].attrs['long_name'] = r'$\mathbf{\mu}_i$'
    return per_state

def get_inter_state(frames: Frames) -> InterState:
    prop: Hashable

    if 'statecomb' in frames:
        #TODO: FIXME: Appropriately remove the 'statecomb' indices and coordinates from the Dataset
        warning(
            "'statecomb' already exists as an index, variable or coordinate"
            " in the dataset, hence it will be removed before recomputation"
        )

    iprops = []
    # TODO: FIXME: check if astate is the correct variable to reference here
    for prop in ['energy', 'nacs', 'astate', 'dip_trans']:
        if prop in frames:
            iprops.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    inter_state = frames[iprops]
    for prop in inter_state:
        if 'state' in inter_state[prop].dims:
            inter_state[prop] = subtract_combinations(
                inter_state[prop], dim='state', labels=True
            )
    inter_state = inter_state.map(keep_norming)

    def state_renamer(lo, hi): 
        if isinstance(lo, int):
            lower_str = f"S_{lo-1}"
        else:
            lower_str = lo
        if isinstance(hi, int):
            higher_str = f"S_{hi-1}"
        else:
            higher_str = hi
        f'${higher_str} - {lower_str}$'

    inter_state = xrhelpers.flatten_midx(
      inter_state,
      'statecomb',
      state_renamer
    )
    if {'energy', 'dip_trans'}.issubset(iprops):
        inter_state = assign_fosc(inter_state)

    inter_state['statecomb'].attrs['long_name'] = "State combinations"
    return inter_state

@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def to_xyz(da: AtXYZ, comment='#') -> str:
    atXYZ = da.transpose('atom', 'direction').values
    atNames = da.atNames.values
    sxyz = np.char.mod('% 23.15f', atXYZ)
    sxyz = np.squeeze(sxyz)
    sxyz = np.hstack((atNames.reshape(-1, 1), sxyz))
    sxyz = np.apply_along_axis(lambda row: ''.join(row), axis=1, arr=sxyz)
    return f'{len(sxyz):>12}\n  {comment}\n' + '\n'.join(sxyz)


@needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
def traj_to_xyz(traj_atXYZ: AtXYZ) -> str:
    atXYZ = traj_atXYZ.transpose(..., 'atom', 'direction').values
    if atXYZ.ndim == 2:
        atXYZ = atXYZ[None, :, :]
    assert len(atXYZ.shape) == 3
    atNames = traj_atXYZ.atNames.values
    sxyz = np.strings.mod('% 13.9f', atXYZ)
    sxyz = atNames[None, :] + sxyz[:, :, 0] + sxyz[:, :, 1] + sxyz[:, :, 2]
    atom_lines = np.broadcast_to([f'{traj_atXYZ.sizes['atom']}'], (sxyz.shape[0], 1))
    if 'time' in traj_atXYZ.coords:
        time_values = np.atleast_1d(traj_atXYZ.coords['time'])
        comment_lines = np.strings.mod('# t=%.2f', time_values)[:, None]
    else:
        comment_lines = np.broadcast_to([''], (sxyz.shape[0], 1))
    return '\n'.join(np.concat([atom_lines, comment_lines, sxyz], 1).ravel())

