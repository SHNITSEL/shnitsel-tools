
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from shnitsel.core.midx import unstack_trajs, stack_trajs


# Function definitions


# Omission strategy
def omit(ds):
    from shnitsel.core.midx import sel_trajs

    return sel_trajs(ds, ds['good_throughout'].all('criterion'))


# Truncation strategy
def cutoffs_to_stacked_mask(ds):
    return ds.coords['time'] < ds['good_upto'].sel(
        trajid_=ds.coords['trajid']
    ).drop_vars('trajid_')


def truncate(ds):
    mask = cutoffs_to_stacked_mask(ds)
    return ds.isel(frame=mask.all('criterion'))


# Transection strategy
def transect(ds, cutoff):
    # TODO: choose cutoff from local maxima of frame counts if not specified
    from shnitsel.core.midx import unstack_trajs

    # If tree:
    if isinstance(ds, xr.DataTree):
        raise NotImplementedError("This function does not yet support DataTree")
    # If stacked:
    elif 'frame' in ds.dims and {'trajid', 'time'}.issubset(ds.coords):
        # stacked = True
        ds = unstack_trajs(ds).fillna(False)
    # Otherwise, had better be unstacked
    elif {'trajid', 'time'}.issubset(ds.dims):
        assert 'is_frame' in ds.coords
    else:
        raise ValueError(
            "The ds argument should be trajectories, " "either stacked or unstacked"
        )
    res = ds.loc[{'time': slice(float(cutoff))}]
    res = res.isel(trajid=res.coords['is_frame'].all('time'))
    res = res.isel(trajid=(res['good_upto'] >= cutoff).all('criterion'))
    return res


def cum_max_quantiles(ds_or_da, quantiles=None):
    from shnitsel.core.midx import unstack_trajs

    if quantiles is None:
        quantiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]

    if hasattr(ds_or_da, 'data_vars'):
        da = ds_or_da['filtranda']
    else:
        da = ds_or_da

    if 'frame' in da.dims:
        da = unstack_trajs(da)

    assert {'trajid', 'time'}.issubset(da.dims)

    da = da.transpose('trajid', 'time', ...).fillna(0)

    cum_max = da.copy(data=np.maximum.accumulate(da.data, axis=1))
    return cum_max.quantile(quantiles, 'trajid')


# Plotting Function!
# TODO: kwargs to modify thresholds
def check_thresholds(ds_or_da, quantiles=None):
    # if 'quantile' in ds_or_da.dims:
    #     assert quantiles is None
    #     quantiles = ds_or_da
    # else:
    #     quantiles = cum_max_quantiles(ds_or_da, quantiles=quantiles)
    if hasattr(ds_or_da, 'data_vars'):
        filtranda = ds_or_da['filtranda']
    else:
        filtranda = ds_or_da
    quantiles = cum_max_quantiles(ds_or_da, quantiles=quantiles)

    fig, axs = plt.subplots(
        quantiles.sizes['criterion'],
        2,
        sharex='col',
        sharey='row',
        layout='constrained',
        width_ratios=[1, 2],
    )
    fig.set_size_inches(6, 2 * quantiles.sizes['criterion'])
    for (title, data), ax in zip(quantiles.groupby('criterion'), axs[:, 1]):
        for qval, qdata in data.groupby('quantile'):
            qdata = qdata.squeeze(['criterion', 'quantile'])
            ax.fill_between(
                qdata.coords['time'], qdata, fc=(0, 0, 0, 0.2), ec=(0, 0, 0, 0)
            )
            ax.text(qdata['time'][-1], qdata[-1], f"{qval*100} %", va='center', c='k')
        if 'thresholds' in data.coords:
            ax.axhline(data.coords['thresholds'].item(), c='r')

    for (title, data), ax in zip(filtranda.groupby('criterion'), axs[:, 0]):
        data = data.squeeze('criterion')
        ax.set_ylabel(title)
        ax.hist(
            data.groupby('trajid').max(),
            density=True,
            cumulative=True,
            orientation='horizontal',
            color='b',
        )
        if 'thresholds' in data.coords:
            ax.axhline(data.coords['thresholds'].item(), c='r')

    axs[-1, 0].set_xlabel('cumulative density\nof per-traj maxima')
    axs[-1, 1].set_xlabel('time / fs')


def cutoffs_to_unstacked_mask(ds, include_total_population=True):
    if 'frame' in ds.dims:
        mask = (ds.coords['time_'] <= ds['good_upto']).rename(
            time_='time', trajid_='trajid'
        )

        if include_total_population:
            is_frame = (
                ds.assign({'is_frame': ('frame', np.ones(frames.sizes['frame']))})[
                    'is_frame'
                ]
                .unstack('frame', fill_value=0)
                .astype(bool)
            )
    else:
        assert {'time', 'trajid'}.issubset(ds.dims)
        mask = ds.coords['time'] < ds['good_upto']

        if include_total_population:
            is_frame = ds.coords['is_frame']

    if include_total_population:
        is_frame = is_frame.expand_dims(
            {'criterion': ['total_population']}
        ).assign_coords(thresholds=('criterion', [np.nan]))
        mask = xr.concat([mask, is_frame], dim='criterion')

    return mask


def validity_populations(ds, intersections=False):
    mask = cutoffs_to_unstacked_mask(ds)
    counts = mask.sum('trajid')
    means = counts.mean('time')
    if intersections:
        counts = mask.sortby(means, ascending=False).cumprod('criterion').sum('trajid')
    else:
        counts = counts.sortby(means, ascending=False)
    fig, axs = plt.subplots(2, 1)
    fig.set_size_inches(6, 8)
    # for label, data in counts.groupby('criterion'):
    for criterion in counts.coords['criterion'].data:
        data = counts.sel(criterion=criterion)
        axs[0].plot(data.coords['time'], data, label=criterion)
        axs[1].plot(data.coords['time'], data * data.coords['time'], label=criterion)
    if intersections:
        order = counts.coords['criterion'].data
        labels = [order[0]] + ['AND ' + x for x in order[1:]]
        axs[0].legend(labels)
    else:
        axs[0].legend()
    axs[0].set_ylabel('# trajectories')
    axs[1].set_ylabel('# frames if transected now')
    axs[1].set_xlabel('time / fs')
    return axs


def energy_filtranda(
    frames,
    *,
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
    from logging import warning
    from shnitsel.core.midx import mdiff
    from shnitsel.units.conversion import convert_energy

    default_thresholds = {
        'etot_drift': etot_drift,
        'etot_step': etot_step,
        'epot_step': epot_step,
        'ekin_step': ekin_step,
        'hop_epot_step': hop_epot_step,
    }

    res = xr.Dataset()
    is_hop = mdiff(frames['astate']) != 0
    e_pot = frames.energy.sel(state=frames.astate).drop_vars('state')
    e_pot.attrs['units'] = frames['energy'].attrs.get('units', 'unknown')
    e_pot = convert_energy(e_pot, to='eV')

    res['epot_step'] = mdiff(e_pot).where(~is_hop, 0)
    res['hop_epot_step'] = mdiff(e_pot).where(is_hop, 0)

    if 'e_kin' in frames.data_vars:
        e_kin = frames['e_kin']
        e_kin.attrs['units'] = frames['e_kin'].attrs.get('units', 'unknown')
        e_kin = convert_energy(e_kin, to='eV')

        e_tot = e_pot + e_kin
        res['etot_drift'] = e_tot.groupby('trajid').map(
            lambda traj: abs(traj - traj.item(0))
        )
        res['ekin_step'] = mdiff(e_kin).where(~is_hop, 0)
        res['etot_step'] = mdiff(e_tot)
    else:
        e_kin = None
        warning("data does not contain kinetic energy variable ('e_kin')")

    da = np.abs(res.to_dataarray('criterion'))
    thresholds = [default_thresholds[x] for x in da.coords['criterion'].data]
    return da.assign_coords(thresholds=('criterion', thresholds))


def last_time_where(mask):
    # mask = mask.transpose('trajid', 'time', ...).astype(bool)
    times = xr.DataArray(np.concat([[-1], mask.coords['time'].data]), dims=('time'))
    good_upto = times[mask.cumprod('time').sum('time')]
    return good_upto


def get_cutoffs(filtranda: xr.DataArray):
    if 'frame' in filtranda.dims and {'trajid', 'time'}.issubset(filtranda.coords):
        stacked = True
        mask = unstack_trajs(filtranda).fillna(False)
    # Otherwise, had better be unstacked
    elif {'trajid', 'time'}.issubset(mask.dims):
        pass  # Format is already as required
    else:
        raise ValueError(
            "The mask argument should be trajectories, " "either stacked or unstacked"
        )
    mask = filtranda < filtranda.coords['thresholds']
    # If tree:

    # tmp_ds = mask.to_dataset('criterion')
    # mask = tmp_ds.assign(total_population=tmp_ds.coords['is_frame']).to_dataarray('criterion')
    mask = mask.astype(bool)

    good_throughout = mask.all('time')
    times = xr.DataArray(np.concat([[-1], mask.coords['time'].data]), dims=('time'))
    good_upto = times[mask.cumprod('time').sum('time')]
    return good_upto.assign_coords(good_throughout)


def assign_cutoffs(ds):
    from shnitsel.core.midx import unstack_trajs

    da = ds['filtranda']
    if isinstance(ds, xr.DataTree):
        raise NotImplementedError("This function does not yet support DataTree")
    # If stacked:
    elif 'frame' in ds.dims and {'trajid', 'time'}.issubset(ds.coords):
        stacked = True
        mask = unstack_trajs(ds).fillna(False)
    # Otherwise, had better be unstacked
    elif {'trajid', 'time'}.issubset(mask.dims):
        pass  # Format is already as required
    else:
        raise ValueError(
            "The mask argument should be trajectories, " "either stacked or unstacked"
        )
    cutoffs = get_cutoffs(da)
    return cutoffs
    # if stacked:
    #     good_throughout = good_throughout.rename(trajid='trajid_')
    #     good_upto = good_upto.rename(trajid='trajid_')
    #     if 'time_' not in ds.coords:
    #         ds = ds.assign_coords(time_=('time_', mask.coords['time'].data))
    # return (
    #     ds
    #     .drop_vars(['filtranda'])
    #     .assign(good_upto=good_upto, good_throughout=good_throughout))


# TODO: plot = True or False; unit conversion
def sanity_check(
    frames,
    cut=False,
    *,
    units='eV',
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
    from numbers import Number

    settings = {
        k: locals()[k]
        for k in ['etot_drift', 'etot_step', 'epot_step', 'ekin_step', 'hop_epot_step']
    }
    frames = frames.assign(filtranda=energy_filtranda(frames, **settings))
    frames = get_cutoffs(frames)
    if not cut:
        return frames
    elif cut == 'truncate':
        return truncate(frames)
    elif cut == 'omit':
        return omit(frames)
    else:
        assert isinstance(cut, Number)
        return transect(frames, cut)
    # if strategy == 'truncate':
    #     return truncate(frames)
    # elif strategy == 'omit':
    #     return omit(frames)
    # elif strategy  == 'transect':
    #     return transect(frames, cutoff)

    # Carolin:
    # if truncate and cutoff is None:
    #     return truncate_fn(frames)
    # elif truncate:
    #     return transect(frames, cutoff)
    # elif not truncate and cutoff is None:
    #     return omit(frames)
    # else:
    #     raise ValueError()
