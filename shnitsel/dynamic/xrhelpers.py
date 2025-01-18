# pyright: basic
import itertools
import os
import xarray as xr
import pandas as pd
from . import postprocess

def midx_combs(values: pd.core.indexes.base.Index|list, name: str|None =None):
    if name is None:
        try:
            # if `values` is a `pandas.core.indexes.base.Index`
            # extract its name
            name = values.name
        except AttributeError:
            raise ValueError("need to specify name if values lack name attribute")

    return xr.Coordinates.from_pandas_multiindex(
      pd.MultiIndex.from_tuples(
        itertools.combinations(values, 2),
        names=['from', 'to']
      ),
      dim=f'{name}comb'
    )


def flatten_midx(ds, idx_name, renamer=None):
    midx = ds.indexes[idx_name]
    to_drop = midx.names + [midx.name]
    fidx = midx.to_flat_index()
    if renamer is not None:
        fidx = [renamer(x, y) for x,y in fidx]
    return ds.drop_vars(to_drop).assign_coords({idx_name: fidx})

def expand_midx(ds, midx_name, level_name, value):
    midx = ds.indexes[midx_name]
    to_drop = [midx.name] + midx.names
    df = midx.to_frame()
    df.insert(0, level_name, [value]*len(midx)) # in place!
    midx = pd.MultiIndex.from_frame(df)
    coords = xr.Coordinates.from_pandas_multiindex(midx, dim=midx_name)
    return (
      ds
      .drop_vars(to_drop)
      .assign_coords(coords)
    )

# TODO! Improve this, then use it to rewrite ts_to_time.
# Only worth doing once you find a second situation in which
# MultiIndex is to be modified.
def modify_midx_level(data, midx_name, coords, replace: set | None):
    """Turns the levels of a MultiIndex into coordinates,
    modifies them as coordinates, then converts back.
    """
    # find MultiIndex name
    # for idx_name, idx in data.indexes.values():
        # if not set(coords).isdisjoint(getattr(idx, 'names', set())):

    # 'frame' in energy_rel.reset_index('frame').coords['ts'].dims

    level_names = list(data.indexes[midx_name].names)
    data = (data
      .reset_index(midx_name)
      .assign_coords(coords)
    )
    if replace is not None:
        new_levels = set(level_names).union(coords.keys()).difference(replace)
        data = (data
          .reset_index(midx_name)
          .set_xindex(list(new_levels))
        )
    return data

def _ts_to_time(data, delta_t=None, swap_index=True):
    default_delta_t = 0.5

    if delta_t is None:
        delta_t = data.attrs.get('delta_t', default_delta_t)

    return modify_midx_level(
      data,
      midx_name='frame',
      coords={'time': data.coords['ts'] * delta_t},
      replace={'ts'} if swap_index else None
    )

def get_frames(path, recode_state=False):
    try:
        frames = xr.open_dataset(path)
    except ValueError as err:
        if not os.path.exists(path):
            raise FileNotFoundError
        else:
            raise err

    tcoord = None
    if 'time' in frames.coords:
        tcoord = 'time'
    elif 'ts' in frames.coords:
        tcoord = 'ts'

    if tcoord is not None:
        frames = frames.set_xindex(['trajid', tcoord]).set_xindex(['from', 'to'])

    if recode_state:
        new_states = [f'$S_{i}$' for i, _ in enumerate(frames['state'])]
        frames = frames.assign_coords({
            'state': new_states,
            'state2': new_states,
            })
    # frames = xrhelpers.modify_midx_level(frames, 'statecomb', {'from':[1,1,2], 'to':[2,3,3]}, )
    frames['energy'] = postprocess.relativize(frames['energy'])
    frames['energy'] = postprocess.convert_energy(frames['energy'], to='eV')
    return frames

def get_frames_time(path):
    return postprocess.ts_to_time(get_frames(path))


def save_frames(frames, path, complevel=9):
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel} for var in frames
    }

    # TODO generalize?
    # netcdf does not support bool
    if 'completed' in frames.data_vars:
        frames = frames.assign(completed=frames.completed.astype('i1'))
    if 'completed' in frames.attrs:
        frames.attrs['completed'] = int(frames.attrs['completed'])
    # or MultiIndex
    frames.reset_index(['frame', 'statecomb']).to_netcdf(
        path, engine='h5netcdf', encoding=encoding
    )
