
import os
import xarray as xr
# TODO: We probably need a fallback version for old and new shnitsel file formats

# def open_frames(path):


def read_shnitsel_file(path: str | os.PathLike):
    """Opens a NetCDF4 file saved by shnitsel-tools, specially interpreting certain attributes.

    Parameters
    ----------
    path
        The path of the file to open.

    Returns
    -------
        An :py:class:`xarray.Dataset` with any MultiIndex restored.

    Raises
    ------
    FileNotFoundError
        If there is is nothing at ``path``, or ``path`` is not a file.
    ValueError (or other exception)
        Raised by the underlying `h5netcdf <https://h5netcdf.org/>`_ engine if the file is corrupted.
    """
    # The error raised for a missing file can be misleading
    try:
        frames = xr.open_dataset(path)
    except ValueError as err:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        else:
            raise err

    # Restore MultiIndexes
    indicator = '_MultiIndex_levels_from_attrs'
    if frames.attrs.get(indicator, False):
        # New way: get level names from attrs
        del frames.attrs[indicator]
        for k, v in frames.attrs.items():
            if k.startswith('_MultiIndex_levels_for_'):
                frames = frames.set_xindex(v)
                del frames.attrs[k]
    else:
        # Old way: hardcoded level names
        tcoord = None
        if 'time' in frames.coords:
            tcoord = 'time'
        elif 'ts' in frames.coords:
            tcoord = 'ts'

        if tcoord is not None:
            frames = frames.set_xindex(['trajid', tcoord])

        if 'from' in frames.coords and 'to' in frames.coords:
            frames = frames.set_xindex(['from', 'to'])

    return frames
