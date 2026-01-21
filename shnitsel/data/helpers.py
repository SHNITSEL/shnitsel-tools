import dataclasses
import logging
import math
from types import UnionType
from typing import Dict, List, Literal, Sequence, TypeVar, get_args
from typing_extensions import TypeForm
import xarray as xr
import numpy as np

from shnitsel.core._api_info import internal
from shnitsel._contracts import needs
from shnitsel.analyze.generic import relativize
from shnitsel.data.multi_indices import sel_trajs
from shnitsel.core.typedefs import DatasetOrArray, Frames
from shnitsel.units.conversion import convert_energy

T = TypeVar("T")

DataType = TypeVar("DataType")
ActualDataType = TypeVar("ActualDataType")


@internal()
def dataclass_from_dict(datatype: type[T], d: List | Dict | T) -> T:
    """Helper function to restore a Dataclass object from its dict representation.

    Mainly used for serialization or storage of data in the DataTree db structure.

    Parameters
    ----------
    datatype : Type[T]
        The dataclass type to restore
    d : List|Dict|T
        The datasource to convert back into the Dataclass instance.

    Raises
    ------
    ValueError
        If value decoding fails during reconstruction
    TypeError
        If some type mismatch occurs between the provided dict and the target type

    Returns
    -------
    T
        A resonstructed `datatype` instance.
    """
    if isinstance(d, list):
        (inner,) = datatype.__args__  # type: ignore
        return [dataclass_from_dict(inner, i) for i in d]  # type: ignore

    try:
        fieldtypes = {f.name: f.type for f in dataclasses.fields(datatype)}  # type: ignore
        return datatype(**{f: dataclass_from_dict(fieldtypes[f], d[f]) for f in d})  # type: ignore
    except Exception as e:
        logging.exception(
            f"Failed to convert object {d} of type {type(d)} to class {datatype}: {e}"
        )
        return d  # type: ignore # Not a dataclass field


def is_assignable_to(
    actual_type: type[ActualDataType] | UnionType,
    expected_dtype: type[DataType] | UnionType | None,
) -> bool:
    """Helper function to check whether a certain type can be considered assignable to
    another type, effectively a more general `is_subclass` function.


    Parameters
    ----------
    actual_type : type[ActualDataType] | TypeForm[ActualDataType]
        The exact type that an object has or that is declared as the type of an object
    expected_dtype : type[DataType] | TypeForm[DataType] | None
        The expected type that we would want to assign to.
        The `actual_type` should be some combination of subclasses of the target types in `expected_dtype`.
        If this parameter is None, the result will be True.

    Returns
    -------
    bool
        True if `actual_type` can be considered some form of subtype to `expected_dtype` or if the
        expected target type is None. Otherwise False.

    """
    if expected_dtype is None:
        return True

    allowed_types: Sequence[type]
    if isinstance(expected_dtype, UnionType):
        allowed_types = list(get_args(expected_dtype))
    else:
        allowed_types = [expected_dtype]

    if isinstance(actual_type, UnionType):
        return all(
            any(issubclass(x, y) for y in allowed_types) for x in get_args(actual_type)
        )
    else:
        return any(issubclass(actual_type, y) for y in allowed_types)


# # TODO: deprecate
# @needs(coords={"ts"})
# def ts_to_time(
#     data: DatasetOrArray,
#     delta_t: float | None = None,
#     old: Literal["drop", "to_var", "keep"] = "drop",
# ) -> DatasetOrArray:
#     assert old in {"drop", "to_var", "keep"}

#     if delta_t is None:
#         if "delta_t" in data:  # could be coord or var
#             # ensure unique
#             arr_delta_t = np.unique(data["delta_t"])
#             assert len(arr_delta_t.shape) == 1
#             if arr_delta_t.shape[0] > 1:
#                 msg = "`delta_t` varies between the trajectories. Please separate the trajectories into groups"
#                 raise ValueError(msg)
#             delta_t = arr_delta_t.item()
#             data = data.drop_vars("delta_t")

#         if "delta_t" in data.attrs:
#             if (
#                 delta_t is not None  # If we already got delta_t from var/coord
#                 and data.attrs["delta_t"] != delta_t
#             ):
#                 msg = "'delta_t' attribute inconsistent with variable/coordinate"
#                 raise ValueError(msg)
#             delta_t = data.attrs["delta_t"]

#         if delta_t is None:  # neither var/coord nor attr
#             msg = "Could not extract `delta_t` from `data`; please pass explicitly"
#             raise ValueError(msg)

#     data = data.reset_index("frame").assign_coords(time=data.coords["ts"] * delta_t)
#     if old in {"drop", "to_var"}:
#         new_levels = list((set(data.indexes["frame"].names) - {"ts"}) | {"time"})
#         data = data.reset_index("frame").set_xindex(new_levels)
#     if old == "drop":
#         data = data.drop_vars("ts")

#     data["time"].attrs.update((dict(units="fs", long_name="$t$", tex_name="t")))
#     data.attrs["delta_t"] = delta_t

#     return data


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
    match to_time, "time" not in ds.coords, "ts" in ds.coords:
        case True, False, _:
            raise ValueError("Timestep coordinate has already been converted to time")
        case True, True, False:
            raise ValueError("No 'ts' coordinate in Dataset")
        case (None, True, True) | (True, True, True):
            ds = ts_to_time(ds)

    match relativize_energy, ds["energy"].min().item() != 0:
        case True, False:
            raise ValueError("Energy is already relativized")
        case (True, True) | (None, True):
            assert "energy" in ds.data_vars
            if relativize_selector is None:
                relativize_selector = {}
            ds = ds.assign({"energy": relativize(ds["energy"], **relativize_selector)})

    match convert_to_eV, ds["energy"].attrs.get("units") != "eV":
        case True, False:
            raise ValueError("Energy is already in eV")
        case (True, True) | (None, True):
            assert "energy" in ds.data_vars
            ds = ds.assign({"energy": convert_energy(ds["energy"], "eV")})

    if convert_e_kin_to_eV and "e_kin" not in ds.data_vars:
        raise ValueError("'frames' object does not have an 'e_kin' variable")
    elif "e_kin" in ds.data_vars:
        match convert_e_kin_to_eV, ds["e_kin"].attrs.get("units") != "eV":
            case True, False:
                raise ValueError("Energy is already in eV")
            case (True, True) | (None, True):
                assert "e_kin" in ds.data_vars
                ds = ds.assign({"e_kin": convert_energy(ds["e_kin"], "eV")})

    return ds


def validate(frames: Frames) -> np.ndarray:
    if "time" in frames.coords:
        tdim = "time"
    elif "ts" in frames.coords:
        tdim = "ts"
    else:
        raise ValueError("Found neither 'time' nor 'ts' coordinate in frames")
    bad_frames = []
    for varname in frames.data_vars.keys():
        # choose appropriate placeholder / bad value for the data_var's dtype
        dtype = frames.dtypes[varname]
        if dtype in {np.dtype("float64"), np.dtype("float32")}:
            mask = np.isnan(frames[varname])
            phname = "`nan`"
        elif dtype in {np.dtype("int32"), np.dtype("int64")}:
            mask = frames[varname] == -1
            phname = "placeholder `-1`"
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
            reddims = set(da.dims) - {"frame"}
            nans = da.sel(frame=mask.any(reddims)).frame
            n = len(nans)
            bfstr = "; ".join(
                [f"trajid={x.trajid.item()} {tdim}={x[tdim].item()}" for x in nans]
            )
            print(f"Variable `{varname}` contains {phname} in {n} frame(s),")
            print(f"    namely: {bfstr}")
            bad_frames += [nans]
        else:
            print(f"Variable `{varname}` does not contain {phname}")

    res: np.ndarray
    if len(bad_frames):
        res = np.unique(xr.concat(bad_frames, dim="frame"))
    else:
        res = np.array([])
    return res


def split_for_saving(frames: Frames, bytes_per_chunk=50e6):
    trajids = frames.get("trajid_", np.unique(frames["trajid"]))
    ntrajs = len(trajids)
    nchunks = math.trunc(frames.nbytes / 50e6)
    logging.debug(f"{nchunks=}")
    indices = np.trunc(np.linspace(0, ntrajs, nchunks + 1)).astype(np.integer)
    logging.debug(f"{indices=}")
    trajidsets = [trajids[a:z].values for a, z in zip(indices[:-1], indices[1:])]
    logging.debug(f"{trajidsets=}")
    return [sel_trajs(frames, trajids[a:z]) for a, z in zip(indices[:-1], indices[1:])]


def save_split(
    frames, path_template, bytes_per_chunk=50e6, complevel=9, ignore_errors=False
):
    from shnitsel.io.shnitsel.write import write_shnitsel_file

    dss = split_for_saving(frames, bytes_per_chunk=bytes_per_chunk)
    for i, ds in enumerate(dss):
        current_path = path_template.format(i)
        try:
            write_shnitsel_file(ds, current_path, complevel=complevel)
        except Exception as e:
            logging.error(f"Exception while saving to {current_path=}")
            if not ignore_errors:
                raise e
