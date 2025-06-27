import numpy as np
import xarray as xr

from .parse.common import __atnum2symbol__
from . import postprocess as P, xrhelpers as xh


def bond_type_to_symbols(e1, e2):
    s1 = __atnum2symbol__[e1]
    s2 = __atnum2symbol__[e2]
    return s1 + s2


def get_bond_types(mol, symbols=True):
    bond_types: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for b in mol.GetBonds():
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()
        indices = (a1.GetIdx(), a2.GetIdx())
        elements = tuple(sorted([a1.GetAtomicNum(), a2.GetAtomicNum()]))
        if elements not in bond_types:
            bond_types[elements] = []
        bond_types[elements].append(indices)
    if symbols:
        return {bond_type_to_symbols(*k): v for k, v in bond_types.items()}
    return bond_types


def get_bond_lengths(atXYZ, bond_types=None, mol=None):
    dists = atXYZ.pipe(P.subtract_combinations, 'atom', labels=True).pipe(P.norm)
    if bond_types is None:
        if mol is None:
            mol = P.default_mol(atXYZ)
        bond_types = get_bond_types(mol, symbols=True)
    return (
        xr.concat(
            [
                dists.sel(atomcomb=bonds).pipe(
                    xh.expand_midx, 'atomcomb', 'bond_type', bond_type
                )
                for bond_type, bonds in bond_types.items()
            ],
            dim='atomcomb',
        )
        .rename({'from': 'atom1', 'to': 'atom2', 'atomcomb': 'bond'})
        .transpose('frame', ...)
    )


def energy_filtranda(frames: xr.Dataset) -> xr.Dataset:
    res = frames[['e_kin']]
    res['e_pot'] = frames.energy.sel(state=frames.astate).drop_vars('state')
    res['e_tot'] = res['e_pot'] + res['e_kin']

    res['etot_drift'] = (
        res['e_tot'].groupby('trajid').map(lambda traj: abs(traj - traj.isel(frame=0)))
    )
    res['ekin_step'] = res['e_kin'].sh.sudi()
    res['epot_step'] = res['e_pot'].sh.sudi()
    res['etot_step'] = res['e_tot'].sh.sudi()
    res['is_hop'] = frames['astate'].sh.sudi() != 0

    return res


def last_time_where(mask):
    mask = mask.unstack('frame', fill_value=False).transpose('trajid', 'time', ...)
    idxs = np.logical_not((~mask.values).cumsum(axis=1)).sum(axis=1)
    times = np.concat([[-1], mask.time.values])
    return mask[:, 0].copy(data=times[idxs]).drop_vars('time').rename('time')


def get_cutoffs(masks_ds):
    ds = masks_ds.map(last_time_where)
    names = list(ds.data_vars)
    ds['original'] = (
        masks_ds.coords['time'].groupby('trajid').max().rename(trajid='trajid_')
    )
    # Put 'original' first, so that min() chooses it in cases of ambiguity
    ds = ds[['original'] + names]
    names = list(ds.data_vars)

    da = ds.to_dataarray('cutoff')
    ds['earliest'] = da.min('cutoff')
    ds['reason'] = da.argmin('cutoff')
    ds.attrs['reasons'] = names
    return ds


def truncate(frames, cutoffs):
    if 'trajid_' not in cutoffs.coords and 'trajid' in cutoffs.coords:
        cutoffs = cutoffs.rename(trajid='trajid_')
    expansion = cutoffs.sel(trajid_=frames.coords['trajid']).drop_vars('trajid_')
    mask = frames['time'] <= expansion
    return frames.sel(frame=mask)


def trajs_where(mask_da): ...