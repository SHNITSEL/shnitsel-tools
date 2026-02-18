from shnitsel import (
    analyze as analyze,
    io as io,
    geo as geo,
    clean as clean,
    filtering as filtering,
    #     vis as vis,
)
# # from shnitsel.data import multi_indices as multi_indices
# # from shnitsel.core.xrhelpers import open_frames as open_frames

from .io import read, write_shnitsel_file, write_ase_db

# # , 'parse', 'open_frames', 'read_trajs', 'read_ase']
# # __all__ = ['io', 'units']
__all__ = [
    'analyze',
    'clean',
    'filtering',
    'io',
    # 'vis',
    # 'data',
    'clean',
    'geo',
    # 'units',
    'read',
    'write_shnitsel_file',
    'write_ase_db',
]


def collapse_display():
    """Collapse or omit verbose representations of Xarray objects"""
    import xarray as xr

    xr.set_options(
        display_expand_coords=False,
        display_expand_data_vars=False,
        display_expand_attrs=False,
        display_expand_data=False,
    )
