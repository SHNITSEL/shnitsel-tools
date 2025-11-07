import math
import itertools
from logging import warning

from typing import Collection, Hashable, Literal, TypeAlias

import numpy as np
import numpy.typing as npt
import xarray as xr

import scipy.stats as st

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401

from .._contracts import needs
from . import xrhelpers

# For backward compatibility while refactoring
from .geom import (
    distance as distance,
    angle as angle,
    dihedral as dihedral,
    normal as normal,
)
from .populations import classical as calc_pops  # noqa: F401
from .ml import pca
from .spectra import (
    assign_fosc,
    get_fosc as get_fosc,
    broaden_gauss as broaden_gauss,
    ds_broaden_gauss as ds_broaden_gauss,
)

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

@needs(dims={'atom'})
def pairwise_dists_pca(atXYZ: AtXYZ, **kwargs) -> xr.DataArray:
    """PCA-reduced pairwise interatomic distances

    Parameters
    ----------
    atXYZ
        A DataArray containing the atomic positions;
        must have a dimension called 'atom'

    Returns
    -------
        A DataArray with the same dimensions as `atXYZ`, except for the 'atom'
        dimension, which is replaced by a dimension 'PC' containing the principal
        components (by default 2)
    """
    res = (
        atXYZ.pipe(subtract_combinations, 'atom')
        .pipe(norm)
        .pipe(pca, 'atomcomb', **kwargs)
    )
    assert not isinstance(res, tuple)  # typing
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


def hop_indices(astates: xr.DataArray) -> xr.DataArray:
    """Find in which frames the active state changes

    Parameters
    ----------
    astates
        A DataArray of state indicators

    Returns
    -------
        A boolean DataArray indicating whether a hop took place
    """
    return sudi(astates) != 0

@needs(coords_or_vars={'atXYZ', 'astate'})
def pca_and_hops(frames: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
    """Get PCA points and info on which of them represent hops

    Parameters
    ----------
    frames
        A Dataset containing 'atXYZ' and 'astate' variables

    Returns
    -------
    pca_res
        The PCA-reduced pairwise interatomic distances
    hops_pca_coords
        `pca_res` filtered by hops, to facilitate marking hops when plotting

    """
    pca_res = pairwise_dists_pca(frames['atXYZ'])
    mask = sudi(frames['astate']) != 0
    hops_pca_coords = pca_res[mask]
    return pca_res, hops_pca_coords

def relativize(da: xr.DataArray, **sel) -> xr.DataArray:
    res = da - da.sel(**sel).min()
    res.attrs = da.attrs
    return res


def setup_frames(
    ds: xr.Dataset,
    *,
    to_time: bool | None = None,
    convert_to_eV: bool | None = None,
    convert_e_kin_to_eV: bool | None = None,
    relativize_energy: bool | None = None,
    relativize_selector=None,
) -> xr.Dataset:
    """Performs several frequent setup tasks.
    Each task can be skipped (by setting the corresponding parameter to False),
    carried out if appropriate (None), or forced in the sense that an error is
    thrown if the task is redundant (True).


    Parameters
    ----------
    ds
        The frames-like xr.Dataset to setup.
    to_time, optional
        Whether to convert a 'ts' (timestep) coordinate to a 'time' coordinate, by default None
    convert_to_eV, optional
        Whether to convert the 'energy' variable to eV, by default None
    convert_e_kin_to_eV, optional
        Whether to convert the 'e_kin' (kinetic energy) variable to eV, by default None
    relativize_energy, optional
        Whether to relativize energies, by default None
    relativize_selector, optional
        This argument is passed to relativize, by default None

    Returns
    -------
        A modified frames-like xr.Dataset

    Raises
    ------
    ValueError
        If a task should be forced (i.e. the corresponding parameter is set to True)
        but cannot be carried out (e.g. because the dataset was already processed previously)
    """
    # TODO: Reconsider how the conversion works here
    match to_time, 'time' not in ds.coords, 'ts' in ds.coords:
        case True, False, _:
            raise ValueError("Timestep coordinate has already been converted to time")
        case True, True, False:
            raise ValueError("No 'ts' coordinate in Dataset")
        case (None, True, True) | (True, True, True):
            ds = ts_to_time(ds)

    match relativize_energy, ds['energy'].min().item() != 0:
        case True, False:
            raise ValueError("Energy is already relativized")
        case (True, True) | (None, True):
            assert 'energy' in ds.data_vars
            if relativize_selector is None:
                relativize_selector = {}
            ds = ds.assign({'energy': relativize(ds['energy'], **relativize_selector)})

    match convert_to_eV, ds['energy'].attrs.get('units') != 'eV':
        case True, False:
            raise ValueError("Energy is already in eV")
        case (True, True) | (None, True):
            assert 'energy' in ds.data_vars
            ds = ds.assign({'energy': convert_energy(ds['energy'], 'eV')})

    if convert_e_kin_to_eV and 'e_kin' not in ds.data_vars:
        raise ValueError("'frames' object does not have an 'e_kin' variable")
    elif 'e_kin' in ds.data_vars:
        match convert_e_kin_to_eV, ds['e_kin'].attrs.get('units') != 'eV':
            case True, False:
                raise ValueError("Energy is already in eV")
            case (True, True) | (None, True):
                assert 'e_kin' in ds.data_vars
                ds = ds.assign({'e_kin': convert_energy(ds['e_kin'], 'eV')})

    return ds



def validate(frames: Frames) -> np.ndarray:
    if 'time' in frames.coords:
        tdim = 'time'
    elif 'ts' in frames.coords:
        tdim = 'ts'
    else:
        raise ValueError("Found neither 'time' nor 'ts' coordinate in frames")
    bad_frames = []
    for varname in frames.data_vars.keys():
        # choose appropriate placeholder / bad value for the data_var's dtype
        dtype = frames.dtypes[varname]
        if dtype in {np.dtype('float64'), np.dtype('float32')}:
            mask = np.isnan(frames[varname])
            phname = '`nan`'
        elif dtype in {np.dtype('int32'), np.dtype('int64')}:
            mask = frames[varname] == -1
            phname = 'placeholder `-1`'
        else:
            print(
                f"Skipping verification of `{varname}` "
                f"as no bad value known for dtype `{dtype}`"
            )

        if mask.all():
            print(
                f"Variable `{varname}` exclusively contains {phname}, "
                "so is effectively missing"
            )
        elif mask.any():
            da = frames[varname]
            reddims = set(da.dims) - {'frame'}
            nans = da.sel(frame=mask.any(reddims)).frame
            n = len(nans)
            bfstr = '; '.join(
                [f"trajid={x.trajid.item()} {tdim}={x[tdim].item()}" for x in nans]
            )
            print(f"Variable `{varname}` contains {phname} in {n} frame(s),")
            print(f"    namely: {bfstr}")
            bad_frames += [nans]
        else:
            print(f"Variable `{varname}` does not contain {phname}")

    res: np.ndarray
    if len(bad_frames):
        res = np.unique(xr.concat(bad_frames, dim='frame'))
    else:
        res = np.array([])
    return res


##############################################
# Functions generally applicable to timeplots:

@needs(coords={'ts'})
def ts_to_time(
    data: xr.Dataset | xr.DataArray,
    delta_t: float | None = None,
    old: Literal['drop', 'to_var', 'keep'] = 'drop',
) -> xr.Dataset | xr.DataArray:
    assert old in {'drop', 'to_var', 'keep'}

    if delta_t is None:
        if 'delta_t' in data:  # could be coord or var
            # ensure unique
            arr_delta_t = np.unique(data['delta_t'])
            assert len(arr_delta_t.shape) == 1
            if arr_delta_t.shape[0] > 1:
                msg = "`delta_t` varies between the trajectories. Please separate the trajectories into groups"
                raise ValueError(msg)
            delta_t = arr_delta_t.item()
            data = data.drop_vars('delta_t')

        if 'delta_t' in data.attrs:
            if (
                delta_t is not None  # If we already got delta_t from var/coord
                and data.attrs['delta_t'] != delta_t
            ):
                msg = "'delta_t' attribute inconsistent with variable/coordinate"
                raise ValueError(msg)
            delta_t = data.attrs['delta_t']

        if delta_t is None:  # neither var/coord nor attr
            msg = "Could not extract `delta_t` from `data`; please pass explicitly"
            raise ValueError(msg)

    data = (data
      .reset_index('frame')
      .assign_coords(time=data.coords['ts'] * delta_t)
    )
    if old in {'drop', 'to_var'}:
        new_levels = list(
            (set(data.indexes['frame'].names) - {'ts'}) | {'time'})
        data = (data
          .reset_index('frame')
          .set_xindex(new_levels)
        )
    if old == 'drop':
        data = data.drop_vars('ts')

    data['time'].attrs.update((dict(units='fs', long_name='$t$', tex_name='t')))
    data.attrs['delta_t'] = delta_t

    return data

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




#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent

def calc_ci(a: npt.NDArray, confidence: float = 0.95) -> npt.NDArray:
    if np.array(a).ndim != 1:
        raise ValueError("This function accepts 1D input only")
    return np.stack(st.t.interval(confidence, len(a)-1, loc=np.mean(a), scale=st.sem(a)))

def ci_agg_last_dim(a, confidence=0.95):
    outer_shape = tuple(a.shape[:-1])
    res = np.full(outer_shape + (3,), np.nan)
    for idxs in np.ndindex(outer_shape):
        res[idxs, :2] = calc_ci(a[idxs], confidence=confidence)
        res[idxs, 2] = np.mean(a[idxs])
    return res

def xr_calc_ci(a: xr.DataArray, dim: DimName, confidence: float = 0.95) -> xr.Dataset:
    res_da: xr.DataArray = xr.apply_ufunc(
        ci_agg_last_dim,
        a,
        kwargs={'confidence': confidence},
        output_core_dims=[['bound']],
        input_core_dims=[[dim]],
    )
    return res_da.assign_coords(  #
        dict(bound=['lower', 'upper', 'mean'])
    ).to_dataset('bound')

@needs(groupable={'time'}, dims={'frame'})
def time_grouped_ci(x: xr.DataArray, confidence: float = 0.9) -> xr.Dataset:
    return (
      x.groupby('time')
      .map(lambda x: xr_calc_ci(x, dim='frame', confidence=confidence)))

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

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index

def set_atom_props(mol, **kws):
    natoms = mol.GetNumAtoms()
    for prop, vals in kws.items():
        if vals is None:
            continue
        elif vals is True:
            vals = range(natoms)
        elif natoms != len(vals):
            raise ValueError(
                f"{len(vals)} values were passed for {prop}, but 'mol' has {natoms} atoms"
            )

        for atom, val in zip(mol.GetAtoms(), vals):
            atom.SetProp(prop, str(val))
    return mol

@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def to_mol(
    atXYZ_frame: xr.DataArray,
    charge: int | None = None,
    covFactor: float = 1.2,
    to2D: bool = True,
    molAtomMapNumber: list | Literal[True] | None = None,
    atomNote: list | Literal[True] | None = None,
    atomLabel: list | Literal[True] | None = None,
) -> rc.Mol:
    """Convert a single frame's geometry to an RDKit Mol object

    Parameters
    ----------
    atXYZ_frame
        The ``xr.DataArray`` object to be converted; must have 'atom' and 'direction' dims,
        must not have 'frame' dim.
    charge
        Charge of the molecule, used by RDKit to determine bond orders; if ``None`` (the default),
        this function will try ``charge=0`` and leave the bond orders undetermined if that causes
        an error; otherwise failure to determine bond order will raise an error.
    covFactor
        Scales the distance at which atoms are considered bonded, by default 1.2
    to2D
        Discard 3D information and generate 2D conformer (useful for displaying), by default True
    molAtomMapNumber
        Set the ``molAtomMapNumber`` properties to values provided in a list,
        or (if ``True`` is passed) set the properties to the respective atom indices
    atomNote
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomNote`` properties
    atomLabel
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomLabel`` properties

    Returns
    -------
        An RDKit Mol object

    Raises
    ------
    ValueError
        If ``charge`` is not ``None`` and bond order determination fails
    """
    mol = rc.rdmolfiles.MolFromXYZBlock(to_xyz(atXYZ_frame))
    rc.rdDetermineBonds.DetermineConnectivity(mol, useVdw=True, covFactor=covFactor)
    try:
        rc.rdDetermineBonds.DetermineBondOrders(mol, charge=(charge or 0))
    except ValueError as err:
        if charge is not None:
            raise err
    if to2D:
        rc.rdDepictor.Compute2DCoords(mol)  # type: ignore
    return set_atom_props(mol, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)


def mol_to_numbered_smiles(mol: rc.Mol) -> str:
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def numbered_smiles_to_mol(smiles: str) -> rc.Mol:
    mol = rc.MolFromSmiles(smiles, sanitize=False)  # sanitizing would strip hydrogens
    map_new_to_old = [-1 for i in range(mol.GetNumAtoms())]
    for atom in mol.GetAtoms():
        # Renumbering with e.g. [3, 2, 0, 1] means atom 3 gets new index 0, not vice-versa!
        map_new_to_old[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    return rc.RenumberAtoms(mol, map_new_to_old)


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def smiles_map(atXYZ_frame, charge=0, covFactor=1.5) -> str:
    mol = to_mol(atXYZ_frame, charge=charge, covFactor=covFactor, to2D=True)
    return mol_to_numbered_smiles(mol)


def default_mol(obj) -> rc.Mol:
    if 'atXYZ' in obj:  # We have a frames Dataset
        atXYZ = obj['atXYZ']
    else:
        atXYZ = obj  # We have an atXYZ DataArray

    if 'smiles_map' in obj.attrs:
        return numbered_smiles_to_mol(obj.attrs['smiles_map'])
    elif 'smiles_map' in atXYZ.attrs:
        return numbered_smiles_to_mol(atXYZ.attrs['smiles_map'])

    try:
        charge = obj.attrs.get('charge', 0)
        return to_mol(atXYZ.isel(frame=0), charge=charge)
    except (KeyError, ValueError):
        raise ValueError(
            "Failed to get default mol, please set a smiles map. "
            "For example, if the compound has charge c and frame i contains a representative geometry, use "
            "frames.attrs['smiles_map'] = frames.atXYZ.isel(frame=i).sh.get_smiles_map(charge=c)"
        )