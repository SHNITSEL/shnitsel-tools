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
from ._version import *

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
    '__version__',
    'show_debug_info'
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

def show_debug_info():
    """Helper function to display version and debugging information of shnitsel (and its dependencies).
    """
    from packaging.version import parse
    from importlib.metadata import version
    from pkgutil import iter_modules
    import sys

    # Get package version
    version_info = parse(__version__)

    # Get python environment version
    python_version = str(sys.version)
    python_version_info = str(sys.version_info)

    # Signal whether we are on a development version
    dev_version_string = "yes" if version_info.dev else "no"

    is_additional_version = version_info.is_devrelease or version_info.is_postrelease or version_info.is_prerelease
    version_suffix_string = ""
    if is_additional_version:
        version_suffix_string = f" (full version: {str(version_info)})"

    # Get environment packages
    packages_str = ""
    for module_info in iter_modules():
        if module_info.ispkg:
            try:
                packages_str += f"{module_info.name} => {version(module_info.name)}\n"
            except: 
                pass
            # packages_str += f"{module_info.name} => [version missing]\n"


    print (f"""
** SHNITSEL TOOLS v{version_info.public}{version_suffix_string} **
Commit: {__commit_id__}
development? {dev_version_string}
Python-version: {python_version}
Python-version_info: {python_version_info}

Loaded Packages:
{packages_str}
           """)