import logging
from types import UnionType
from typing import TYPE_CHECKING, Any, Iterable, Sequence, TypeVar, overload
from typing_extensions import TypeForm

import numpy as np

from shnitsel.core._api_info import API, internal
from shnitsel.data.dataset_containers import Trajectory, Frames, wrap_dataset
import xarray as xr

from shnitsel.data.dataset_containers.shared import ShnitselDataset

if TYPE_CHECKING:
    from shnitsel.data.tree import ShnitselDB

# TODO: FIXME: Set units on delta_t and t_max when converted into a variable

_coordinate_meta_keys = ["trajid", "delta_t", "max_ts", "t_max", "completed", "nsteps"]


class InconsistentAttributeError(ValueError):
    pass


class MultipleCompoundsError(ValueError):
    pass


class MissingValue:
    """Sentinel value for ``tree_to_frames``."""


@internal()
def _check_matching_dimensions(
    datasets: Iterable[xr.Dataset | Trajectory | Frames | xr.DataArray],
    excluded_dimensions: set[str] = set(),
    limited_dimensions: set[str] | None = None,
) -> bool:
    """Function to check whether all/certain dimensions are equally sized.

    Excluded dimensions can be provided as a set of strings.

    Parameters
    ----------
    datasets : Iterable[xr.Dataset]
        The series of datasets to be checked for equal dimensions
    excluded_dimensions : set[str], optional
        The set of dimension names to be excluded from the comparison. Defaults to set().
    limited_dimensions : set[str], optional
        Optionally set a list of dimensions to which the analysis should be limited.

    Returns
    -------
    bool
        True if all non-excluded (possibly limited) dimensions match in size.  False otherwise.
    """
    # TODO: FIXME: Should we check that the values are also the same?

    res_matching = True
    matching_dims = {}
    distinct_dims = []
    is_first = True

    for ds in datasets:
        for dim in ds.dims:
            if str(dim) in excluded_dimensions:
                # Do not bother with excluded dimensions
                continue

            if limited_dimensions is not None and str(dim) not in limited_dimensions:
                # Skip if we are not in the set list of limited_dimensions
                continue

            if is_first:
                matching_dims[str(dim)] = ds.sizes[dim]
            else:
                if (
                    str(dim) not in matching_dims
                    or matching_dims[str(dim)] != ds.sizes[dim]
                ):
                    res_matching = False
                    distinct_dims.append(str(dim))
        is_first = False

    logging.info(f"Found discrepancies in the following dimensions: {distinct_dims}")

    return res_matching


@internal()
def _compare_dicts_of_values(
    curr_root_a: Any, curr_root_b: Any, base_key: list[str] = []
) -> tuple[list[list[str]] | None, list[list[str]] | None]:
    """Compare two dicts and return the lists of matching and non-matching recursive keys.

    Parameters
    ----------
    curr_root_a : Any
        Root of the first tree
    curr_root_b : Any
        Root of the second tree
    base_key : list[str]
        The current key associated with the root. Starts with [] for the initial call.

    Returns
    -------
    tuple[list[list[str]] | None, list[list[str]] | None]
        A tuple, where the first list is the list of chains of keys of all matching sub-trees,
        the second entry is the same but for identifying distinct sub-trees.
        If a matching key points to a sub-tree, the entire sub-tree is identical.
    """
    matching_keys = []
    non_matching_keys = []
    if curr_root_a == curr_root_b:
        # This subtree matches
        return ([base_key], None)
    else:
        if isinstance(curr_root_a, dict) and isinstance(curr_root_b, dict):
            # We need to recurse further
            keys_a = set(curr_root_a.keys())
            keys_b = set(curr_root_a.keys())
            delta_keys = keys_a.symmetric_difference(keys_b)
            shared_keys = keys_a.intersection(keys_b)

            for key in delta_keys:
                non_matching_keys.append(base_key + [key])

            for key in shared_keys:
                new_base = base_key + [key]

                if key not in curr_root_a or key not in curr_root_b:
                    non_matching_keys.append(new_base)
                    continue

                res_matching, res_non_matching = _compare_dicts_of_values(
                    curr_root_a[key], curr_root_b[key], new_base
                )

                if res_matching is not None:
                    matching_keys.extend(res_matching)
                if res_non_matching is not None:
                    non_matching_keys.extend(res_non_matching)

            return (
                None if len(matching_keys) == 0 else matching_keys,
                None if len(non_matching_keys) == 0 else non_matching_keys,
            )
        else:
            # This subtree does not match and we do not need to recurse further
            return (None, [base_key])


@internal()
def _check_matching_var_meta(
    datasets: Sequence[xr.Dataset | Trajectory | Frames | xr.DataArray],
) -> bool:
    """Function to check if all of the variables have matching metadata.

    We do not want to merge trajectories with different metadata on variables.

    TODO: Allow for variables being denoted that we do not care for.

    Parameters
    ----------
    datasets : Sequence[xr.Dataset | Trajectory | Frames]
        The trajectories to compare the variable metadata for.

    Returns
    -------
    bool
        True if the metadata matches on all trajectories, False otherwise
    """
    collected_meta = []

    shared_vars = None

    for ds in datasets:
        ds_vars = ds.coords if isinstance(ds, xr.DataArray) else ds.variables
        ds_meta = {}
        this_vars = set(ds_vars.keys())
        if shared_vars is None:
            shared_vars = this_vars
        else:
            shared_vars = this_vars.intersection(shared_vars)

        for var_name in ds_vars:
            var_attr = ds[var_name].attrs.copy()
            ds_meta[var_name] = var_attr
        collected_meta.append(ds_meta)

    if shared_vars is None:
        return True

    # TODO: FIXME: This should probably fail if variables are not present on all datasets.

    for i in range(len(datasets) - 1):
        for var in shared_vars:
            _matching, distinct_keys = _compare_dicts_of_values(
                collected_meta[i][var], collected_meta[i + 1][var]
            )

            if distinct_keys is not None and len(distinct_keys) > 0:
                return False

    return True


@internal()
def _merge_traj_metadata(
    datasets: Sequence[xr.Dataset | Trajectory | Frames | xr.DataArray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Function to gather metadate from a set of trajectories.

    Used to combine trajectories into one aggregate Dataset.

    Parameters
    ----------
    datasets : Sequence[xr.Dataset | Trajectory | Frames]
        The sequence of trajctories for which metadata should be collected

    Returns
    -------
    tuple[dict[str, Any], dict[str, np.ndarray]]
            The resulting meta information shared across all trajectories (first),
            and then the distinct meta information (second) in a key -> Array_of_values fashion.
    """
    num_datasets = len(datasets)
    shared_meta = {}
    distinct_meta = {}

    if num_datasets == 0:
        return shared_meta, distinct_meta

    traj_meta_distinct_defaults = {
        "__mol": np.full((num_datasets,), None, dtype="O"),
        "trajid": np.full((num_datasets,), -1, dtype="i4"),
        "delta_t": np.full((num_datasets,), np.nan, dtype="f8"),
        "max_ts": np.full((num_datasets,), -1, dtype="i4"),
        "t_max": np.full((num_datasets,), np.nan, dtype="f8"),
        "completed": np.full((num_datasets,), False, dtype="?"),
        "nsteps": np.full((num_datasets,), -1, dtype="i4"),
    }

    # Assert the existence of a trajectory id for each trajectory.
    all_keys = set()
    all_keys.add("trajid")

    for ds in datasets:
        for x in ds.attrs.keys():
            x_str = str(x)
            if not x_str.startswith("__"):
                # ignore private attrs
                all_keys.add(str(x))

    all_keys.intersection_update([str(k) for k in traj_meta_distinct_defaults.keys()])

    all_meta = {}
    for key in all_keys:
        kept_array = None
        if key in traj_meta_distinct_defaults:
            kept_array = traj_meta_distinct_defaults[key]
        else:
            kept_array = np.full((num_datasets,), None, dtype=object)

        for i, ds in enumerate(datasets):
            if key in ds.attrs:
                kept_array[i] = ds.attrs[key]
            else:
                if isinstance(ds, ShnitselDataset):
                    if key == 'trajid':
                        kept_array[i] = ds.trajectory_id
                    elif key == 'delta_t':
                        kept_array[i] = ds.delta_t
                    elif key == 't_max':
                        kept_array[i] = ds.t_max
                elif isinstance(ds, xr.DataArray):
                    # Special treatment for data arrays necessary
                    if key == 'trajid':
                        kept_array[i] = (
                            ds.coords['trajectory'].item()
                            if 'trajectory' in ds.coords
                            else (
                                ds.coords['atrajectory'].min()
                                if 'atrajectory' in ds.coords
                                else ds.attrs.get(
                                    "trajid",
                                    ds.attrs.get(
                                        "trajectory", ds.attrs.get("trajectory_id", -1)
                                    ),
                                )
                            )
                        )
                    elif key == 'delta_t':
                        kept_array[i] = (
                            ds.coords['delta_t'].item()
                            if 'delta_t' in ds.coords
                            else -1
                        )
                    elif key == 't_max':
                        kept_array[i] = (
                            ds.coords['t_max'].item() if 't_max' in ds.coords else -1
                        )

        all_meta[key] = kept_array

    keep_distinct = ["__mol", "trajid", "delta_t", "max_ts", "t_max", "completed"]

    for key in all_keys:
        if key in keep_distinct:
            # We treat some specific values different
            distinct_meta[key] = all_meta[key]
        else:
            try:
                set_of_vals = set(all_meta[key])

                # If there are distinct meta values, we assign the values all to the distinct set. Otherwise, we only keep the one as shared.
                if len(set_of_vals) > 1:
                    distinct_meta[key] = all_meta[key]
                else:
                    shared_meta[key] = set_of_vals.pop()
            except TypeError:
                distinct_meta[key] = all_meta[key]

    # Add missing trajectory ids and reassign duplicate ids
    used_trajectory_ids = set()
    next_candidate_id = 0

    for i in range(num_datasets):
        if (
            distinct_meta["trajid"][i] < 0
            or distinct_meta["trajid"][i] is None
            or distinct_meta["trajid"][i] in used_trajectory_ids
        ):
            while next_candidate_id in used_trajectory_ids:
                next_candidate_id += 1
            distinct_meta["trajid"][i] = next_candidate_id

    return shared_meta, distinct_meta


DataType = TypeVar("DataType")


@overload
@API()
def concat_trajs(
    datasets: Sequence[xr.DataArray],
    dtype: type[DataType] | UnionType | None = None,
) -> xr.DataArray: ...


@overload
@API()
def concat_trajs(
    datasets: Sequence[Trajectory | Frames | xr.Dataset],
    dtype: type[DataType] | UnionType | None = None,
) -> xr.Dataset: ...


@API()
def concat_trajs(
    datasets: Sequence[Trajectory | Frames | xr.Dataset] | Sequence[xr.DataArray],
    dtype: type[DataType] | UnionType | None = None,
) -> xr.Dataset | xr.DataArray:
    """Function to concatenate multiple trajectories along their `time` dimension.

    Will create one continuous time dimension like an extended trajectory.
    The concatenated dimension will be renamed `frame` consisting of a `time` and a `atrajectory` component
    where the latter denotes the active trajectory.

    Additionally, a dimension `trajectory` with accompanying trajectory ids as metadata and
    to index the remaining collected trajectory metadata will be introduced.

    For a sequence of data arrays, we will just try and concatenate the arrays.

    Parameters
    ----------
    datasets : Iterable[Trajectory | Frames | xr.Dataset] | Sequence[xr.DataArray]
        Datasets representing the individual trajectories or a sequence of arrays to concatenate.
    dtype :  type[DataType] | UnionType | None
        Type hint for the data to be included in the resulting container type.

    Raises
    ------
    ValueError
        Raised if there is conflicting input dimensions.
    ValueError
        Raised if there is conflicting input variable meta data.
    ValueError
        Raised if there is conflicting global input attributes that are relevant to the merging process.
    ValueError
        Raised if there are no trajectories provided to this function.

    Returns
    -------
    xr.Dataset
        The combined and extended trajectory with a new leading `frame` dimension
    """

    all_traj = all(isinstance(x, (Trajectory, Frames, xr.Dataset)) for x in datasets)
    all_da = all(isinstance(x, (xr.DataArray)) for x in datasets)

    if not all_traj and not all_da:
        logging.error(
            "Attempted to concatenate mixed trajectory and non-trajectory data"
        )
        raise ValueError(
            "Concatenation is only possible for pure Trajectory/xarray.Dataset or pure xarray.DataArray input data. The provided data contained unsupported or mixed types."
        )

    if all_traj:
        wrapped_ds = [wrap_dataset(x, (Trajectory | Frames)) for x in datasets]

        datasets_pure = list(
            x.as_frames if isinstance(x, Trajectory) else x for x in wrapped_ds
        )

        if len(datasets_pure) == 0:
            raise ValueError("No trajectories were provided.")

        is_multi_trajectory = len(datasets_pure) > 1

        # Check that we do not have pre-existing multi-trajectory datasets
        for ds in datasets_pure:
            if ds.is_multi_trajectory:
                is_multi_trajectory = True
                # TODO: FIXME: This should actually not be a problem. Look into this
                logging.error(
                    "Multi-trajectory dataset provided to concat() function. Aborting."
                )
                raise ValueError(
                    "Multi-trajectory dataset provided to concat() function."
                )

            # TODO: FIXME: Deal with multi-trajectory merges
            # if (
            #     "trajid" in ds.coords
            #     or "trajid_" in ds.coords
            #     or "frame" in ds.coords
            #     or 'atrajectory' in ds.coords
            #     or 'trajectory' in ds.coords
            #     or 'trajectory' in ds.dims
            # ):
            #     dsid = (
            #         ds.coords["trajid"]
            #         if "trajid" in ds.coords
            #         else ds.attrs["trajid"]
            #         if "trajid" in ds.attrs
            #         else "unknown"
            #     )
            #     raise ValueError(
            #         f"trajectory with existing `trajid`={dsid} (or `trajid_`), 'atrajectory', 'trajectory' or existing `frame` coordinates provided to `concat`. Indicates prior merge of multiple trajetories. Cannot proceed. Please only provide trajectories without these coordinates"
            #     )

        # Check that all dimensions match. May want to check the values match as well?
        if not _check_matching_dimensions(datasets_pure, {"frame"}):
            message = "Dimensions of the provided data vary."
            logging.warning(
                f"{message} Merge result may be inconsistent. Please ensure you only merge consistent trajectories."
            )
            # TODO: Do we want to merge anyway?
            raise ValueError(f"{message} Will not merge.")

        # All units should be converted to same unit
        if not _check_matching_var_meta(datasets_pure):
            # TODO: FIXME: Add message info which variable did not match.
            message = (
                "Variable meta attributes vary between different tajectories. "
                "This indicates inconsistencies like distinct units between trajectories. "
                "Please ensure consistency between datasets before merging."
            )
            logging.warning(f"{message} Merge result may be inconsistent.")
            # TODO: Do we want to merge anyway?
            raise ValueError(f"{message} Will not merge.")

        # trajid set by merge_traj_metadata
        consistent_metadata, distinct_metadata = _merge_traj_metadata(datasets_pure)

        # To keep trajid as a part of the multi-index distinct from
        datasets_amended = [
            # Expansion applied in `as_frames` property call.
            # ds.expand_dims(atrajectory=[distinct_metadata["trajid"][i]]).stack(
            #     frame=["atrajectory", "time"]
            # )
            ds.dataset
            for i, ds in enumerate(datasets_pure)
        ]

        # TODO: Check if the order of datasets stays the same. Otherwise distinct attributes may not be appropriately sorted.
        frames = xr.concat(
            datasets_amended, dim="frame", coords="different", combine_attrs="override"
        )

        # DataArrays without a `trajectory` dimension cannot have these coords assigned
        # Introduce new trajid dimension
        frames = frames.assign_coords(
            trajectory=(
                "trajectory",
                distinct_metadata["trajid"],
                {"description": "id of the original trajectory before concatenation."},
            )
        )

        # Add remaining trajectory-metadata
        # First the ones that may end up as coordinates
        frames = frames.assign_coords(
            {
                k: (
                    "trajectory",
                    v,
                    {"description:": f"Attribute {k} merged in concatenation"},
                )
                for k, v in distinct_metadata.items()
                if k != "trajid" and str(k) in _coordinate_meta_keys
            }
        )
    elif all_da:
        # Concatenate data arrays

        def da_to_frame(da, default_id):
            if 'frame' in da.dims:
                return da

            if 'time' not in da.dims:
                raise ValueError(
                    f"The data array {da=} did not have sufficient information for concatenation. Missing `time` dimension."
                )

            if 'atrajectory' not in da.coords:
                trajectory_id = da.attrs.get(
                    "trajid",
                    da.attrs.get("trajectory_id", da.attrs.get("trajectory", None)),
                )
                if trajectory_id is None:
                    for source in [
                        'trajid',
                        'atrajectory',
                        'trajectory',
                        'trajectory_id',
                    ]:
                        if source in da.coords:
                            trajectory_id = da.coords[source].item()
                if trajectory_id is None:
                    trajectory_id = default_id
                    # raise ValueError(
                    #     f"Data Array {da=} did not contain any trajectory id information. Cannot concatenate."
                    # )

                da = da.expand_dims(atrajectory=[trajectory_id])

            return da.stack(frame=["atrajectory", "time"])

        framed_da: list[xr.DataArray] = [
            da_to_frame(da, i) for i, da in enumerate(datasets)
        ]

        # Check that all dimensions match. May want to check the values match as well?
        if not _check_matching_dimensions(framed_da, {"frame"}):
            message = "Dimensions of the provided data vary."
            logging.warning(
                f"{message} Merge result may be inconsistent. Please ensure you only merge consistent trajectories."
            )
            # TODO: Do we want to merge anyway?
            raise ValueError(f"{message} Will not merge.")

        # All units should be converted to same unit
        if not _check_matching_var_meta(framed_da):
            # TODO: FIXME: Add message info which variable did not match.
            message = (
                "Variable meta attributes vary between different tajectories. "
                "This indicates inconsistencies like distinct units between trajectories. "
                "Please ensure consistency between datasets before merging."
            )
            logging.warning(f"{message} Merge result may be inconsistent.")
            # TODO: Do we want to merge anyway?
            raise ValueError(f"{message} Will not merge.")

        # trajid set by merge_traj_metadata
        consistent_metadata, distinct_metadata = _merge_traj_metadata(framed_da)

        # TODO: Check if the order of datasets stays the same. Otherwise distinct attributes may not be appropriately sorted.
        # TODO: FIXME: This has some FutureWarnings associated with it by xarray. Keep an eye on this.
        frames = xr.concat(framed_da, dim="frame", coords='different', compat="equals")
        is_multi_trajectory = True
    else:
        raise RuntimeError("Something went wrong in data concatenation")

    # Set merged metadata
    frames.attrs.update(consistent_metadata)

    # Previous update
    # frames = frames.assign_coords(trajid_=traj_meta["trajid"])
    # frames = frames.assign(
    #    delta_t=("trajid", traj_meta["delta_t"]),
    #    max_ts=("trajid", traj_meta["max_ts"]),
    #    completed=("trajid", traj_meta["completed"]),
    #    nsteps=("trajid", traj_meta["nsteps"]),
    # )

    # Then all remaining metadata
    frames.attrs.update(
        {
            k: v
            for k, v in distinct_metadata.items()
            if k != "trajid" and str(k) not in _coordinate_meta_keys
        }
    )

    # Envelop in the wrapper proxy
    # if not isinstance(frames, Trajectory):
    #     frames = Trajectory(frames)

    # if TYPE_CHECKING:
    #     assert isinstance(frames, Trajectory)

    frames.attrs["is_multi_trajectory"] = is_multi_trajectory

    return frames


@API()
def db_from_data(
    datasets: Sequence[DataType] | DataType,
    dtype: type[DataType] | UnionType | None = None,
) -> "ShnitselDB[DataType]":  # noqa: F821
    """Function to merge multiple trajectories of the same molecule into a single ShnitselDB instance.

    Parameters
    ----------
    datasets : Sequence[DataType] | DataType
        The individual loaded data points, e.g. trajectories or a single data point/trajectory to turn into a tree.
    dtype :  type[DataType] | UnionType | None
        Type hint for the data to be included in the resulting tree.

    Returns
    -------
    ShnitselDB[DataType]
        The resulting ShnitselDB structure with a ShnitselDBRoot, CompoundGroup and DataGroup layers.
    """
    from shnitsel.data.tree import ShnitselDB, complete_shnitsel_tree

    if isinstance(datasets, Sequence):
        if len(datasets) == 0:
            return ShnitselDB[DataType](dtype=dtype)

        # Collect trajectories, check if trajectories match and build databases
        datasets_list = list(datasets)
        # TODO: FIXME: This should probably require more throrough testing?
        # And potentially derivation of different compounds?

        if isinstance(datasets_list[0], (Trajectory, Frames)):
            if not _check_matching_dimensions(
                [x.dataset for x in datasets_list],  # type: ignore # entries should be all trajectories or frames here or this all breaks anyway
                limited_dimensions=set("atom"),
            ):
                raise ValueError(
                    "Could not merge trajectories into one ShnitselDB, because compound `unknown` would contain distinct compounds. "
                    "Please only load one type of compound at a time."
                )

        return complete_shnitsel_tree(datasets_list)
    else:
        # We only need to wrap a single trajectory
        return complete_shnitsel_tree(datasets)


@API()
def layer_trajs(
    datasets: Sequence[xr.Dataset | Trajectory | Frames],
    dtype: type[DataType] | UnionType | None = None,
) -> xr.Dataset:
    """Function to combine trajectories into one Dataset by creating a new dimension 'trajectory' and indexing the different trajectories along that.

    Will create one new trajectory dimension.

    Parameters
    ----------
    datasets : Sequence[xr.Dataset | Trajectory]
        Datasets representing the individual trajectories
    dtype :  type[DataType] | UnionType | None
        Type hint for the data to be included in the resulting container type.

    Raises:
    ValueError
        Raised if there is conflicting input meta data.
    ValueError
        Raised if there are no trajectories provided to this function or if there are non-trajectories provided to this function.


    Returns
    -------
    xr.Dataset
        The combined and extended trajectory with a new leading `trajectory` dimension to differentiate the trajectory data.
    """

    datasets_converted = list(
        x.dataset if isinstance(x, ShnitselDataset) else x for x in datasets
    )

    if len(datasets_converted) == 0:
        raise ValueError("No trajectories were provided.")

    if not _check_matching_dimensions(datasets_converted, {'time', 'frame'}):
        message = "Dimensions of the provided datasets are not consistent."
        logging.warning(
            f"{message} Merge result may be inconsistent. Please ensure you only merge consistent trajectories."
        )
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    # All units should be converted to same unit
    if not _check_matching_var_meta(datasets_converted):
        # TODO: FIXME: Add message info which variable did not match.
        message = (
            "Variable meta attributes vary between different tajectories. "
            "This indicates inconsistencies like distinct units between trajectories. "
            "Please ensure consistency between datasets before merging."
        )
        logging.warning(f"{message} Merge result may be inconsistent.")
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    consistent_metadata, distinct_metadata = _merge_traj_metadata(datasets_converted)

    trajids = distinct_metadata["trajid"]

    datasets = [
        ds.expand_dims(trajectory=[id]) for ds, id in zip(datasets_converted, trajids)
    ]

    # trajids = pd.Index(meta["trajid"], name="trajid")
    # coords_trajids = xr.Coordinates(indexes={"trajid": trajids})
    # breakpoint()
    layers = xr.concat(
        datasets, dim="trajectory", combine_attrs="drop_conflicts", join="outer"
    )

    # layers = layers.assign_coords(trajid=trajids)

    # del meta["trajid"]
    # layers = layers.assign(
    #    {k: xr.DataArray(v, dims=["trajid"])
    #     for k, v in meta.items() if k != "trajid"}
    # )
    layers.attrs.update(consistent_metadata)

    # Add remaining trajectory-metadata
    layers = layers.assign_coords(
        {
            k: (
                "trajectory",
                v,
                {"description:": f"Attribute {k} merged in concatenation"},
            )
            for k, v in distinct_metadata.items()
            if k != "trajid" and str(k) in _coordinate_meta_keys
        }
    )

    # Then all remaining metadata
    layers.attrs.update(
        {
            k: v
            for k, v in distinct_metadata.items()
            if k != "trajid" and str(k) not in _coordinate_meta_keys
        }
    )

    layers.attrs["is_multi_trajectory"] = True

    if not isinstance(layers, xr.Dataset):
        layers = xr.Dataset(layers)

    if TYPE_CHECKING:
        assert isinstance(layers, xr.Dataset)

    return layers
