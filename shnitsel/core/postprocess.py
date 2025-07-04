import math
import itertools
from logging import warning

from typing import Collection, Hashable, Literal, TypeAlias

import numpy as np
import xarray as xr

import scipy.stats as st
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
import rdkit.Chem as rc

from . import xrhelpers

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
    """Calculate the 2-norm of a DataArray, reducing the dimension with dame `dim`

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
        A DataArray with dimension `dim` reduced
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
        xrmat = xr.DataArray(
            data=mat, coords={f'{dim}comb': midx(da, dim), dim: da.indexes[dim]}
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

def pca_for_plot(diffnorms):
    """Legacy method to calculate 2-component PCA

    Parameters
    ----------
    diffnorms
        Data to be transformed, usually norms of pairwise differences

    Returns
    -------
    data
        The transformed data
    pca_n2_scaled
        The trained PCA object produced by scikit-learn
    """
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.decomposition import PCA
    
    
    scaled = MinMaxScaler()\
        .fit_transform(diffnorms)
    pca_n2_scaled = PCA(n_components=2)
    pca_n2_scaled.fit(scaled)
    return pca_n2_scaled.transform(scaled), pca_n2_scaled

def pca(
    da: xr.DataArray, dim: str, n_components: int = 2, return_pca_object=False
) -> tuple[xr.DataArray, PCA] | xr.DataArray:
    """xarray-oriented wrapper around scikit-learn's PCA

    Parameters
    ----------
    da
        A DataArray with at least a dimension with a name matching `dim`
    dim
        The name of the dimension to reduce
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False

    Returns
    -------
        A DataArray with the same dimensions as `da`, except for the dimension
        indicated by `dim`, which is replaced by a dimension `PC` of size `n_components`
    """
    scaled = xr.apply_ufunc(
      MinMaxScaler().fit_transform,
      da.transpose(..., dim)
    )
    
    pca_object = PCA(n_components=n_components)
    pca_object.fit(scaled)
    pca_res: xr.DataArray = xr.apply_ufunc(
        pca_object.transform,
        scaled,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )

    if return_pca_object:
        return (pca_res, pca_object)
    else:
        return pca_res

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

def _sudi_groupby(da):
    """Successive differences"""
    da = da.transpose('frame', ...)
    return da.groupby('trajid').map(
        lambda traj: np.diff(
            traj, axis=0, prepend=np.array(traj[0], ndmin=traj.values.ndim)
        )
    )

def sudi(da: xr.DataArray) -> xr.DataArray:
    """Successive differences"""
    da = da.transpose('frame', ...)
    res = np.diff(da, axis=0, prepend=np.array(da[0], ndmin=da.values.ndim))
    # Don't compare the last timestep of one trajectory to the first timestep of the next:
    res[da.time == 0] = 0
    return da.copy(data=res)


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
    axidx_frame = astates.get_axis_num("frame")
    assert isinstance(axidx_frame, int)
    conseq_diffs = np.diff(astates, axis=axidx_frame, prepend=0)
    return astates.copy(data=conseq_diffs) != 0


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
    mask = hop_indices(frames['astate'])
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


def convert(
    da: xr.DataArray, to: str, quantity: str, conversions: dict
) -> xr.DataArray:
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
        res: xr.DataArray = da * dividend / divisor
    res.attrs.update({'units': to})
    return res

class Converter:
    def __init__(self, quantity, conversions):
        self.quantity = quantity
        self.conversions = conversions
        self.targets = list(self.conversions.keys())

    def __call__(self, da: xr.DataArray, to: str) -> xr.DataArray:
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
            res: xr.DataArray = da * dividend / divisor
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

convert_length = Converter(
    'length',
    dict(
        pm=52.91772105,
        angstrom=0.5291772105,
        bohr=1,
    ),
)

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

def ts_to_time(
    data: xr.Dataset | xr.DataArray,
    delta_t: float | None = None,
    old: Literal['drop', 'to_var', 'keep'] = 'drop',
):
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


# def aggregate_trajs(frames):
#     gb = frames.groupby('trajid')
#     return gb.mean(), gb.stddev() 

# def state_diffs(prpt):
#     statecomb = xrhelpers.get_statecombs


def _get_fosc(energy, dip_trans):
    return 2 / 3 * energy * dip_trans**2


def assign_fosc(ds: xr.Dataset) -> xr.Dataset:
    da = _get_fosc(convert_energy(ds['energy'], to='hartree'), ds['dip_trans'])
    da.name = 'fosc'
    da.attrs['long_name'] = r"$f_{\mathrm{osc}}$"
    return ds.assign(fosc=da)

def broaden_gauss(
    E: xr.DataArray,
    fosc: xr.DataArray,
    agg_dim: DimName = 'frame',
    *,
    width: float = 0.5,
    nsamples: int = 1000,
    xmax: float | None = None,
) -> xr.DataArray:
    r"""
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
    res: xr.DataArray = (g(Espace - E) * fosc).mean(dim=agg_dim)
    res.name = 'fosc'
    res.attrs = fosc.attrs
    for cname, coord in res.coords.items():
        if cname in fosc.coords:
            coord.attrs = fosc.coords[cname].attrs
    return res.assign_coords({'energy': Espace})

def ds_broaden_gauss(
    ds: xr.Dataset, width: float = 0.5, nsamples: int = 1000, xmax: float | None = None
) -> xr.DataArray:
    return broaden_gauss(
        ds['energy'], ds['fosc'], width=width, nsamples=nsamples, xmax=None
    )


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
    iprops = []
    for prop in ['energy', 'nacs', 'astate', 'dip_trans']:
        if prop in frames:
            iprops.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    inter_state = frames[iprops]
    for prop in inter_state:
        if 'state' in inter_state[prop].dims:
            inter_state[prop] = subtract_combinations(inter_state[prop], dim='state')

    inter_state = inter_state.map(keep_norming)
    inter_state = xrhelpers.flatten_midx(
      inter_state,
      'statecomb',
      lambda lo, hi: f'$S_{hi-1} - S_{lo-1}$'
    )
    if {'energy', 'dip_trans'}.issubset(iprops):
        inter_state = assign_fosc(inter_state)

    inter_state['statecomb'].attrs['long_name'] = "States"
    return inter_state

def calc_pops(frames: Frames) -> xr.DataArray:
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
    nstates = frames.sizes['state']
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
    return (pops / pops.sum('state')).assign_coords(state=frames['state'])


#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent

def calc_ci(a, confidence=0.95):
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


def time_grouped_ci(x: xr.DataArray, confidence: float = 0.9) -> xr.Dataset:
    return (
      x.groupby('time')
      .map(lambda x: xr_calc_ci(x, dim='frame', confidence=confidence)))

def to_xyz(da: AtXYZ, comment='#'):
    atXYZ = da.values
    atNames = da.atNames.values
    sxyz = np.char.mod('% 23.15f', atXYZ)
    sxyz = np.squeeze(sxyz)
    sxyz = np.hstack((atNames.reshape(-1, 1), sxyz))
    sxyz = np.apply_along_axis(lambda row: ''.join(row), axis=1, arr=sxyz)
    return f'{len(sxyz):>12}\n  {comment}\n' + '\n'.join(sxyz)


def traj_to_xyz(traj_atXYZ: AtXYZ):
    return '\n'.join(
        to_xyz(t_atXYZ, comment=f"# t={t}") for t, t_atXYZ in traj_atXYZ.groupby('time')
    )


######################################################
# Functions relating to calculation of dihedral angles
def dnorm(a): return norm(a, dim='direction')
def dcross(a, b): return xr.cross(a, b, dim='direction')
def ddot(a, b): return xr.dot(a, b, dim='direction')
def angle_(a, b): return np.arccos(ddot(a, b) / (dnorm(a) * dnorm(b)))
def normal(a, b, c): return dcross(a-b, c-b)

def dihedral_(a, b, c, d):
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    return angle_(abc, bcd)

def full_dihedral_(a, b, c, d):
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    sign = np.sign(ddot(dcross(abc, bcd), (c - b)))
    return sign * angle_(abc, bcd)


def dihedral(
    atXYZ: AtXYZ,
    i: int,
    j: int,
    k: int,
    l: int,
    *,
    deg: bool = False,
    full: bool = False,
) -> xr.DataArray:
    a = atXYZ.isel(atom=i)    
    b = atXYZ.isel(atom=j)
    c = atXYZ.isel(atom=k)
    d = atXYZ.isel(atom=l)
    result: xr.DataArray = full_dihedral_(a, b, c, d) if full else dihedral_(a, b, c, d)
    if deg:
        result = result * 180 / np.pi
    result.name = 'dihedral'
    result.attrs['long_name'] = r"$\varphi_{%d,%d,%d,%d}$" % (i, j, k, l)
    return result


def angle(atXYZ: AtXYZ, i: int, j: int, k: int, *, deg: bool = False) -> xr.DataArray:
    a = atXYZ.isel(atom=i)    
    b = atXYZ.isel(atom=j)
    c = atXYZ.isel(atom=k)
    ab = a-b
    cb = c-b
    result: xr.DataArray = angle_(ab, cb)
    if deg:
        result = result * 180 / np.pi
    result.name = 'angle'
    result.attrs['long_name'] = r"$\theta_{%d,%d,%d}$" % (i, j, k)
    return result


def distance(atXYZ: AtXYZ, i: int, j: int) -> xr.DataArray:
    a = atXYZ.isel(atom=i)    
    b = atXYZ.isel(atom=j)
    result: xr.DataArray = dnorm(a - b)
    result.name = 'distance'
    result.attrs['long_name'] = r"$\|\mathbf{r}_{%d,%d}\|$" % (i, j)
    return result


###############################################
# Functions to investigate hops in a trajectory
# Note: some of these functions represent statecombs
# using complex numbers, because MultiIndex was
# getting awkward

def trajs_with_hops(astates: Astates) -> list[Hashable]:
    """Example usage: `trajs_with_hops(frames['astate'])`
    """
    return [
      trajid for trajid, traj in astates.groupby('trajid')
      if len(np.unique(traj)) > 1]

def get_hop_types(astates: Astates) -> dict[int, tuple[int, int]]:
    """Example usage:
    """
    pairs = np.c_[astates[:-1], astates[1:]]
    hop_types = {}
    for i, (s1, s2) in enumerate(pairs):
        if s1 != s2:
            hop_types[i] = (s1, s2)
    return hop_types

def pick_statecombs(
    da: xr.DataArray,
    statecombs: xr.DataArray,
    frames: Frames,
    framedim: DimName = 'frame',
) -> xr.DataArray:
    assert len(statecombs) == len(frames)
    if 'statecomb' not in da.sizes:
        # no picking to do
        return da.isel({framedim: frames}).copy()
    # translate statecombs labels to indices
    picks = {
        'statecomb': da.indexes['statecomb'].get_indexer(statecombs),
        framedim: frames}
    # but not frames, as these should be indices already

    tmp = [
        range(size) if (x := picks.get(dim)) is None else x
        for dim, size in da.sizes.items()
    ]
    indexer = tuple(tmp)

    coords = da.isel({framedim: frames}).coords.copy()
    del(coords['statecomb'])

    return xr.DataArray(da.values[indexer], coords=coords)

def find_traj_hops(traj: xr.Dataset) -> xr.Dataset:
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

    return traj.map(
        pick_statecombs, statecombs=statecombs, frames=frames, framedim=framedim
    ).assign(statecomb=xr.DataArray(statecombs, dims=[framedim]))

def find_hops(frames: Frames) -> Frames:
    mask = frames['trajid'].isin(trajs_with_hops(frames['astate']))
    return (
      frames.sel(frame=mask)
      .assign_coords(statecomb=[1+2j, 1+3j, 2+3j]) # TODO generalize... and/or sanitize inputs!
      .groupby('trajid').map(find_traj_hops))

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index

def to_mol(atXYZ_frame, charge=None, covFactor=1.5, to2D=True):
    mol = rc.rdmolfiles.MolFromXYZBlock(to_xyz(atXYZ_frame))
    rc.rdDetermineBonds.DetermineConnectivity(mol, useVdw=True, covFactor=covFactor)
    try:
        rc.rdDetermineBonds.DetermineBondOrders(mol, charge=(charge or 0))
    except ValueError as err:
        if charge is not None:
            raise err
    if to2D:
        rc.rdDepictor.Compute2DCoords(mol)  # type: ignore
    return mol


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

def default_mol(obj):
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