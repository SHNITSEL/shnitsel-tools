import math
import itertools
from logging import warning

import numpy as np
import xarray as xr

import scipy.stats as st
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
import rdkit.Chem as rc

from . import xrhelpers

_var_delta_t_msg = "`delta_t` varies between the trajectories. Please separate the trajectories into groups"


def norm(da, dim='direction', keep_attrs=None):
    return xr.apply_ufunc(
        np.linalg.norm,
        da,
        input_core_dims=[[dim]],
        on_missing_core_dim='copy',
        kwargs={"axis": -1},
        keep_attrs=keep_attrs
    )

def subtract_combinations(da, dim, labels=False):
    def midx(da, dim):
        return xrhelpers.midx_combs(da.indexes[dim])[f'{dim}comb']

    if dim not in da.dims:
        raise ValueError(f"'{dim}' is not a dimension of the DataArray")
    
    n = da.sizes[dim]

    mat = np.zeros((math.comb(n, 2), n))
    combs = itertools.combinations(range(n), 2)

    # After matrix multiplication, index r of output vector has value c2 - c1
    for r, (c1, c2) in enumerate(combs):
        mat[r, c1] = -1
        mat[r, c2] = 1

    if labels:
        mat = xr.DataArray(
            data=mat,
            coords={
              f'{dim}comb': midx(da, dim),
              dim: da.indexes[dim]
            }
        )
    else:
        mat = xr.DataArray(
            data=mat,
            dims=[f'{dim}comb', dim],
        )

    newdims = list(da.dims)
    newdims[newdims.index(dim)] = f'{dim}comb'

    res = (mat @ da).transpose(*newdims)
    res.attrs = da.attrs
    res.attrs['deltaed'] = set(res.attrs.get('deltaed', [])).union({dim})
    return res

def pca_for_plot(diffnorms):
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.decomposition import PCA
    
    
    scaled = MinMaxScaler()\
        .fit_transform(diffnorms)
    pca_n2_scaled = PCA(n_components=2)
    pca_n2_scaled.fit(scaled)
    return pca_n2_scaled.transform(scaled), pca_n2_scaled

def pca(
  da: xr.DataArray,
  dim: str,
  n_components: int=2,
  return_pca_object=False
) -> tuple[xr.DataArray, PCA] | xr.DataArray:
    scaled = xr.apply_ufunc(
      MinMaxScaler().fit_transform,
      da.transpose(..., dim)
    )
    
    pca_object = PCA(n_components=n_components)
    pca_object.fit(scaled)
    pca_res = xr.apply_ufunc(
      pca_object.transform,
      scaled,
      input_core_dims=[[dim]],
      output_core_dims=[['PC']]
    )

    if return_pca_object:
        return (pca_res, pca_object)
    else:
        return pca_res

def pairwise_dists_pca(atXYZ, **kwargs):
    return (atXYZ
      .pipe(subtract_combinations, 'atom')
      .pipe(norm)
      .pipe(pca, 'atomcomb', **kwargs))

def hop_indices(astates):
    axidx_frame: int = astates.get_axis_num("frame")
    conseq_diffs = np.diff(astates, axis=axidx_frame, prepend=0)
    return astates.copy(data=conseq_diffs) != 0     

def pca_and_hops(frames: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
    pca_res = pairwise_dists_pca(frames['atXYZ'])
    mask = hop_indices(frames['astate'])
    hops_pca_coords = pca_res[mask]
    return pca_res, hops_pca_coords

def relativize(da: xr.DataArray, **sel):
    res = da - da.sel(**sel).min()
    res.attrs = da.attrs
    return res

def convert(da: xr.DataArray, to: str, quantity: str, conversions: dict):
    try:
        from_ = da.attrs['units']
    except AttributeError:
        raise TypeError("da should be a DataArray with a da.attr attribute.")
    except KeyError:
        raise KeyError("The 'units' attribute of the DataArray must be set.")

    try:
        divisor = conversions[from_]
    except KeyError:
        targets = list(conversions.keys())
        raise ValueError(f"Can't convert {quantity} from {from_!r}, only from: {targets}")

    try:
        dividend = conversions[to]
    except KeyError:
        targets = list(conversions.keys())
        raise ValueError(f"Can't convert {quantity} to {to!r}, only to: {targets}")

    with xr.set_options(keep_attrs=True):
        res = da * dividend / divisor
    res.attrs.update({'units': to})
    return res

class Converter:
    def __init__(self, quantity, conversions):
        self.quantity = quantity
        self.conversions = conversions
        self.targets = list(self.conversions.keys())


    def __call__(self, da: xr.DataArray, to: str):
        try:
            from_ = da.attrs['units']
        except AttributeError:
            raise TypeError("da should be a DataArray with a da.attr attribute.")
        except KeyError:
            raise KeyError("The 'units' attribute of the DataArray must be set.")

        try:
            divisor = self.conversions[from_]
        except KeyError:
            raise ValueError(f"Can't convert {self.quantity} from {from_!r}, only from: {self.targets}")
    
        try:
            dividend = self.conversions[to]
        except KeyError:
            raise ValueError(f"Can't convert {self.quantity} to {to!r}, only to: {self.targets}")
    
        with xr.set_options(keep_attrs=True):
            res = da * dividend / divisor
        res.attrs.update({'units': to})
        return res

convert_energy = Converter('energy', dict(
  hartree=1.0,
  au=1.0,
  eV=27.211386245988,
  keV=0.027211386245988
))

convert_dipoles = Converter('dipoles', dict(
  au=1.0,
  debye=1/0.3934303
))

# def convert_energy(da: xr.DataArray, to: str):
#     conversions = dict(
#         hartree=1.0,
#         eV=27.211386245988,
#         keV=0.027211386245988
#     )
#     return convert(da, to, quantity='energy', conversions=conversions)

def changes(da):
    # TODO
    diffs = da.copy(data=np.diff(
      da.coords['astates'],
      axis=da.coords['astates'].get_axis_num("frame"),
      prepend=0) 
    )
    return diffs != 0

# def get_hops(ds, )

def validate(frames) -> np.ndarray:
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

    if len(bad_frames):
        return np.unique(xr.concat(bad_frames, dim='frame'))
    else:
        return np.array([])


##############################################
# Functions generally applicable to timeplots:

def ts_to_time(data, delta_t=None, old='drop'):
    assert old in {'drop', 'to_var', 'keep'}

    if delta_t is None:
        if 'delta_t' in data:  # could be coord or var
            # ensure unique
            if (data['delta_t'] == data['delta_t'][0]).all():
                delta_t = data['delta_t'].item(0)
                if data.attrs.get('delta_t', delta_t) != delta_t:
                    raise ValueError(
                        "'delta_t' attribute inconsistent with variable/coordinate"
                    )
                data = data.drop_vars('delta_t')
                data.attrs['delta_t'] = delta_t
            else:
                raise ValueError(_var_delta_t_msg)
        elif 'delta_t' in data.attrs:
            delta_t = data.attrs['delta_t']
        else:
            raise ValueError(
                "Could not extract `delta_t` from `data`; please pass explicitly"
            )

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

    return data

def keep_norming(da, exclude=None):
    if exclude is None:
        exclude = {'state', 'statecomb', 'frame'}
    for dim in set(da.dims) - exclude:
        da = norm(da, dim, keep_attrs=True)
        da.attrs['norm_order'] = 2
    return da


# def aggregate_trajs(frames):
#     gb = frames.groupby('trajid')
#     return gb.mean(), gb.stddev() 

# def state_diffs(prpt):
#     statecomb = xrhelpers.get_statecombs


def _get_fosc(energy, dip_trans):
    return 2 / 3 * energy * dip_trans**2


def assign_fosc(ds):
    da = _get_fosc(convert_energy(ds['energy'], to='hartree'), ds['dip_trans'])
    da.name = 'fosc'
    da.attrs['long_name'] = r"$f_{\mathrm{osc}}$"
    return ds.assign(fosc=da)

def broaden_gauss(E, fosc, agg_dim='frame', *, width=0.5, nsamples=1000, xmax=None):
    """
    Parameters
    ----------
    E
        values used for the x-axis, presumably $E_i$
    fosc
        values used for the y-axis, presumably $f_\mathrm{osc}$
    agg_dim, optional
        dimension along which to aggregate the many Gaussian distributions,
        by default 'frame'
    width, optional
        the width (i.e. 2 standard deviations) of the Gaussian distributions
        used, by default 0.001
    nsamples, optional
        number of evenly spaced x-values over which to sample the distribution,
        by default 1000
    xmax, optional
        the maximum x-value, by default 3 standard deviations
        beyond the pre-broadened maximum
    """

    stdev = width / 2

    def g(x):
        nonlocal stdev
        return 1 / (np.sqrt(2 * np.pi) * stdev) * np.exp(-(x**2) / (2 * stdev**2))

    if xmax is None:
        # broadening could visibly overshoot the former maximum by 3 standard deviations
        xmax = E.max().item() * (1 + 1.5 * width)
    xs = np.linspace(0, xmax, num=nsamples)
    Espace = xr.DataArray(
        xs,
        dims=['energy'],
        attrs=E.attrs)
    res = (g(Espace - E) * fosc).mean(dim=agg_dim)
    res.name = 'fosc'
    res.attrs = fosc.attrs
    for cname, coord in res.coords.items():
        if cname in fosc.coords:
            coord.attrs = fosc.coords[cname].attrs
    return res.assign_coords({'energy': Espace})

def ds_broaden_gauss(ds, width=0.001, nsamples=1000, xmax=None):
    return broaden_gauss(
        ds['energy'], ds['fosc'], width=width, nsamples=nsamples, xmax=None
    )


def get_per_state(frames):
    props_per = {'energy', 'forces', 'dip_perm'}.intersection(frames.keys())
    per_state = frames[props_per].map(keep_norming, keep_attrs=False)
    per_state['forces'] = per_state['forces'].where(per_state['forces'] != 0)

    per_state['energy'].attrs['long_name'] = r'$E$'
    per_state['forces'].attrs['long_name'] = r'$\mathbf{F}$'
    if 'dip_perm' in per_state:
        per_state['dip_perm'].attrs['long_name'] = r'$\mathbf{\mu}_i$'
    return per_state

def get_inter_state(frames):
    iprops = []
    for prop in ['energy', 'nacs', 'astate', 'dip_trans']:
        if prop in frames:
            iprops.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    inter_state = frames[iprops]
    inter_state['energy'] = subtract_combinations(inter_state['energy'], dim='state')

    inter_state = inter_state.apply(keep_norming)
    inter_state = xrhelpers.flatten_midx(
      inter_state,
      'statecomb',
      lambda lo, hi: f'$S_{hi-1} - S_{lo-1}$'
    )
    inter_state['statecomb'].attrs['long_name'] = "States"
    return inter_state

def calc_pops(frames):
    """Fast way to calculate populations
    Requires states ids to be small integers
    """
    data = frames['astate']
    if -1 in frames['astate']:
        warning(
            "`frames['astate']` contains the placeholder value `-1`, "
            "indicating missing state information.  "
            "The frames in question will be excluded from the "
            "population count altogether."
        )
        data = data.sel(frame=(data != -1))
    nstates=frames.sizes['state']
    # zero_or_one = int(frames.coords['state'].min())
    zero_or_one = 1  # TODO: For now, assume lowest state is 1
    assert zero_or_one in {0,1}
    pops = data.groupby('time').map(
        lambda group: xr.apply_ufunc(
            lambda values: np.bincount(values, minlength=nstates + zero_or_one)[
                zero_or_one:
            ],
            group,
            input_core_dims=[['frame']],
            output_core_dims=[['state']],
        )
    )
    return pops / pops.sum('state')

#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent

def calc_ci(a, confidence=0.95):
    if np.array(a).ndim != 1:
        raise ValueError("This function accepts 1D input only")
    return np.stack(st.t.interval(confidence, len(a)-1, loc=np.mean(a), scale=st.sem(a)))

def ci_agg_last_dim(a, confidence=0.95):
    outer_shape = a.shape[:-1]
    res = np.full(outer_shape + (3,), np.nan)
    for idxs in np.ndindex(outer_shape):
        res[idxs, :2] = calc_ci(a[idxs], confidence=confidence)
        res[idxs, 2] = np.mean(a[idxs])
    return res

def xr_calc_ci(a, dim, confidence=0.95):
    return xr.apply_ufunc(
        ci_agg_last_dim,
        a,
        kwargs={'confidence': confidence},
        output_core_dims=[['bound']],
        input_core_dims=[[dim]]
    ).assign_coords(dict(bound=['lower', 'upper', 'mean'])
    ).to_dataset('bound')

def time_grouped_ci(x, confidence=0.9):
    return (
      x.groupby('time')
      .map(lambda x: xr_calc_ci(x, dim='frame', confidence=confidence)))

def to_xyz(da, comment='#'):
    atXYZ = da.values
    atNames = da.atNames.values
    sxyz = np.char.mod('%s', atXYZ)
    sxyz = np.squeeze(sxyz)
    sxyz = np.hstack((atNames.reshape(-1, 1), sxyz))
    sxyz = np.apply_along_axis(lambda row: '\t'.join(row), axis=1, arr=sxyz)
    return f'{len(sxyz):>12}\n  {comment}\n' + '\n'.join(sxyz)


def traj_to_xyz(traj_atXYZ):
    return '\n'.join(
        to_xyz(t_atXYZ, comment=f"# t={t}") for t, t_atXYZ in traj_atXYZ.groupby('time')
    )


######################################################
# Functions relating to calculation of dihedral angles
def dnorm(a): return norm(a, dim='direction')
def dcross(a, b): return xr.cross(a, b, dim='direction')
def ddot(a, b): return xr.dot(a, b, dim='direction')
def angle(a, b): return np.arccos(ddot(a, b) / (dnorm(a) * dnorm(b)))
def normal(a, b, c): return dcross(a-b, c-b)

def dihedral_(a, b, c, d):
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    return angle(abc, bcd)

def dihedral(atXYZ, i, j, k, l):
    a = atXYZ.isel(atom=i)    
    b = atXYZ.isel(atom=j)
    c = atXYZ.isel(atom=k)
    d = atXYZ.isel(atom=l)
    data = dihedral_(a, b, c, d)
    return data

def distance(atXYZ, i, j):
    a = atXYZ.isel(atom=i)    
    b = atXYZ.isel(atom=j)
    data = dnorm(a - b)
    return data

###############################################
# Functions to investigate hops in a trajectory
# Note: some of these functions represent statecombs
# using complex numbers, because MultiIndex was
# getting awkward

def trajs_with_hops(astates):
    """Example usage: `trajs_with_hops(frames['astate'])`
    """
    return [
      trajid for trajid, traj in astates.groupby('trajid')
      if len(np.unique(traj)) > 1]

def get_hop_types(astates):
    """Example usage:
    """
    pairs = np.c_[astates[:-1], astates[1:]]
    hop_types = {}
    for i, (s1, s2) in enumerate(pairs):
        if s1 != s2:
            hop_types[i] = (s1, s2)
    return hop_types

def pick_statecombs(da, statecombs, frames, framedim='frame'):
    assert len(statecombs) == len(frames)
    if 'statecomb' not in da.sizes:
        # no picking to do
        return da.isel({framedim: frames}).copy()
    # translate statecombs labels to indices
    picks = {
        'statecomb': da.indexes['statecomb'].get_indexer(statecombs),
        framedim: frames}
    # but not frames, as these should be indices already

    indexer = [
        range(size) if (x:=picks.get(dim)) is None else x
        for dim, size in da.sizes.items()]
    
    coords = da.isel({framedim: frames}).coords.copy()
    del(coords['statecomb'])

    return xr.DataArray(da.values[*indexer], coords=coords)

def find_traj_hops(traj):
    def check(s): return s if s in traj.sizes else False
    framedim = check('frame') or check('time') or 'ts'

    hops = get_hop_types(traj['astate'])
    if len(hops) == 0:
        return (
          traj.isel({framedim: [0]})
          .sel({'statecomb': []}, drop=True)
          # .assign(statecomb=xr.DataArray([np.nan], dims=[framedim]))
        )

    frames, statecombs = [], []
    for idx, h in hops.items():
        frames += [idx, idx+1]
        statecombs += [min(h)+max(h)*1j]*2

    return (
      traj.apply(pick_statecombs, statecombs=statecombs, frames=frames, framedim=framedim)
      .assign(statecomb=xr.DataArray(statecombs, dims=[framedim])))

def find_hops(frames):
    mask = frames['trajid'].isin(trajs_with_hops(frames['astate']))
    return (
      frames.sel(frame=mask)
      .assign_coords(statecomb=[1+2j, 1+3j, 2+3j]) # TODO generalize... and/or sanitize inputs!
      .groupby('trajid').map(find_traj_hops))

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index


def mol_to_numbered_smiles(mol):
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def numbered_smiles_to_mol(smiles):
    mol = rc.MolFromSmiles(smiles, sanitize=False)  # sanitizing would strip hydrogens
    map_new_to_old = [-1 for i in range(mol.GetNumAtoms())]
    for atom in mol.GetAtoms():
        # Renumbering with e.g. [3, 2, 0, 1] means atom 3 gets new index 0, not vice-versa!
        map_new_to_old[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    return rc.RenumberAtoms(mol, map_new_to_old)
