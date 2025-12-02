import logging
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Set, Tuple

import numpy as np

from shnitsel.data.shnitsel_db_format import ShnitselDB, build_shnitsel_db
from shnitsel.data.trajectory_format import Trajectory
import xarray as xr

# TODO: FIXME: Set units on delta_t and t_max when converted into a variable

def _check_matching_dimensions(
    datasets: Iterable[Trajectory],
    excluded_dimensions: Set[str] = set(),
    limited_dimensions: Set[str] | None = None,
) -> bool:
    """Function to check whether all/certain dimensions are equally sized.

    Excluded dimensions can be provided as a set of strings.

    Args:
        datasets (Iterable[Trajectory]): The series of datasets to be checked for equal dimensions
        excluded_dimensions (Set[str], optional): The set of dimension names to be excluded from the comparison. Defaults to set().
        limited_dimensions (Set[str], optional): Optionally set a list of dimensions to which the analysis should be limited.

    Returns:
        bool: True if all non-excluded (possibly limited) dimensions match in size.  False otherwise.
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


def _compare_dicts_of_values(
    curr_root_a: Any, curr_root_b: Any, base_key: List[str] = []
) -> Tuple[List[List[str]] | None, List[List[str]] | None]:
    """Compare two dicts and return the lists of matching and non-matching recursive keys.

    Args:
        curr_root_a (Any): Root of the first tree
        curr_root_b (Any): Root of the second tree
        base_key (List[str]): The current key associated with the root. Starts with [] for the initial call.

    Returns:
        Tuple[List[List[str]]|None, List[List[str]]|None]: A tuple, where the first list is the list of chains of keys of all matching sub-trees,
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


def _check_matching_var_meta(
    datasets: List[Trajectory],
) -> bool:
    """Function to check if all of the variables have matching metadata.

    We do not want to merge trajectories with different metadata on variables.

    TODO: Allow for variables being denoted that we do not care for.

    Args:
        datasets (List[Trajectory]): The trajectories to compare the variable metadata for.

    Returns:
        bool: True if the metadata matches on all trajectories, False otherwise
    """
    collected_meta = []

    shared_vars = None

    for ds in datasets:
        ds_meta = {}
        this_vars = set(ds.variables.keys())
        if shared_vars is None:
            shared_vars = this_vars
        else:
            shared_vars = this_vars.intersection(shared_vars)

        for var_name in ds.variables:
            var_attr = ds[var_name].attrs.copy()
            ds_meta[var_name] = var_attr
        collected_meta.append(ds_meta)

    for i in range(len(datasets) - 1):
        for var in shared_vars:
            _matching, distinct_keys = _compare_dicts_of_values(
                collected_meta[i][var], collected_meta[i + 1][var]
            )

            if distinct_keys is not None and len(distinct_keys) > 0:
                return False

    return True


def _merge_traj_metadata(
    datasets: List[Trajectory],
) -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    """Function to gather metadate from a set of trajectories.

    Used to combine trajectories into one aggregate Dataset.

    Args:
        datasets (Iterable[Trajectory]): The sequence of trajctories for which metadata should be collected

    Returns:
        Tuple[Dict[str,Any],Dict[str,np.ndarray]]: The resulting meta information shared across all trajectories (first),
                and then the distinct meta information (second) in a key -> Array_of_values fashion.
    """
    num_datasets = len(datasets)
    shared_meta = {}
    distinct_meta = {}

    if num_datasets == 0:
        return shared_meta, distinct_meta

    traj_meta_distinct_defaults = {
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

        all_meta[key] = kept_array

    keep_distinct = ["trajid", "delta_t", "max_ts", "t_max", "completed"]

    for key in all_keys:
        if key in keep_distinct:
            # We treat some specific values different
            distinct_meta[key] = all_meta[key]
        else:
            set_of_vals = set(all_meta[key])

            # If there are distinct meta values, we assign the values all to the distinct set. Otherwise, we only keep the one as shared.
            if len(set_of_vals) > 1:
                distinct_meta[key] = all_meta[key]
            else:
                shared_meta[key] = set_of_vals.pop()

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


def concat_trajs(datasets: Iterable[Trajectory]) -> Trajectory:
    """Function to concatenate multiple trajectories along their `time` dimension.

    Will create one continuous time dimension like an extended trajectory. The concatenated dimension will be renamed `frame`

    Args:
        datasets (Iterable[Trajectory]): Datasets representing the individual trajectories

    Raises:
        ValueError: Raised if there is conflicting input dimensions.
        ValueError: Raised if there is conflicting input variable meta data.
        ValueError: Raised if there is conflicting global input attributes that are relevant to the merging process.
        ValueError: Raised if there are no trajectories provided to this function.

    Returns:
        Trajectory: The combined and extended trajectory with a new leading `frame` dimension
    """

    datasets = list(datasets)

    if len(datasets) == 0:
        raise ValueError("No trajectories were provided.")

    # Check that we do not have pre-existing multi-trajectory datasets
    for ds in datasets:
        if "is_multi_trajectory" in ds.attrs:
            logging.error(
                "Multi-trajectory dataset provide to concat() function. Aborting."
            )
            raise ValueError("Multi-trajectory dataset provided to concat() function.")

        if "trajid" in ds.coords or "trajid_" in ds.coords or "frame" in ds.coords:
            dsid = (
                ds.coords["trajid"]
                if "trajid" in ds.coords
                else ds.attrs["trajid"]
                if "trajid" in ds.attrs
                else "unknown"
            )
            raise ValueError(
                f"trajectory with existing `trajid`={dsid} (or `trajid_`) or existing `frame` coordinates provided to `concat`. Indicates prior merge of multiple trajetories. Cannot proceed. Please only provide trajectories without these coordinates"
            )

    # Check that all dimensions match. May want to check the values match as well?
    if not _check_matching_dimensions(datasets, set(["time"])):
        message = "Dimensions of the provided data vary."
        logging.warning(
            f"{message} Merge result may be inconsistent. Please ensure you only merge consistent trajectories."
        )
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    # All units should be converted to same unit
    if not _check_matching_var_meta(datasets):
        message = (
            "Variable meta attributes vary between different tajectories. "
            "This indicates inconsistencies like distinct units between trajectories. "
            "Please ensure consistency between datasets before merging."
        )
        logging.warning(f"{message} Merge result may be inconsistent.")
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    # trajid set by merge_traj_metadata
    consistent_metadata, distinct_metadata = _merge_traj_metadata(datasets)

    # To keep trajid as a part of the multi-index distinct from
    datasets = [
        ds.expand_dims(trajid_=[distinct_metadata["trajid"][i]]).stack(
            frame=["trajid_", "time"]
        )
        for i, ds in enumerate(datasets)
    ]

    # TODO: Check if the order of datasets stays the same. Otherwise distinct attributes may not be appropriately sorted.
    frames = xr.concat(datasets, dim="frame", combine_attrs="drop_conflicts")

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

    # Introduce new trajid dimension
    frames = frames.assign_coords(
        trajid=(
            "trajid",
            distinct_metadata["trajid"],
            {"description": "id of the original trajectory before concatenation."},
        )
    )

    # Add remaining trajectory-metadata
    frames = frames.assign_coords(
        {
            k: ("trajid", v, {"description:": f"Attribute {k} merged in concatenation"})
            for k, v in distinct_metadata
            if k != "trajid"
        }
    )

    # Envelop in the wrapper proxy
    if not isinstance(frames, Trajectory):
        frames = Trajectory(frames)

    if TYPE_CHECKING:
        assert isinstance(frames, Trajectory)

    frames.attrs["is_multi_trajectory"] = True

    return frames


def db_from_trajs(datasets: Iterable[Trajectory] | Trajectory) -> ShnitselDB:  # noqa: F821
    """Function to merge multiple trajectories of the same molecule into a single ShnitselDB instance.

    Args:
        datasets (Iterable[Trajectory]): The individual loaded trajectories.

    Returns:
        ShnitselDB: The resulting ShnitselDB structure with a ShnitselDBRoot, CompoundGroup and TrajectoryData layers.
    """
    if not isinstance(datasets, Trajectory):
        # Collect trajectories, check if trajectories match and build databases
        datasets_list = list(datasets)
        if not _check_matching_dimensions(
            datasets_list, limited_dimensions=set("atom")
        ):
            raise ValueError(
                "Could not merge datasets into one ShnitselDB, because compound `unknown` would contain distinct compounds. "
                "Please only load one type of compound at a time."
            )

        return build_shnitsel_db(datasets_list)
    else:
        # We only need to wrap a single trajectory
        return build_shnitsel_db(datasets)


def layer_trajs(datasets: Iterable[Trajectory]) -> Trajectory:
    """Function to combine trajectories into one Dataset by creating a new dimension 'trajid' and indexing the different trajectories along that.

    Will create one new trajid dimension.

    Args:
        datasets (Iterable[xr.Dataset]): Datasets representing the individual trajectories

    Raises:
        ValueError: Raised if there is conflicting input meta data.
        ValueError: Raised if there are no trajectories provided to this function.


    Returns:
        xr.Dataset: The combined and extended trajectory with a new leading `trajid` dimension
    """

    datasets = list(datasets)

    if len(datasets) == 0:
        raise ValueError("No trajectories were provided.")

    if not _check_matching_dimensions(datasets, set(["time"])):
        message = "Dimensions of the provided data vary."
        logging.warning(
            f"{message} Merge result may be inconsistent. Please ensure you only merge consistent trajectories."
        )
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    # All units should be converted to same unit
    if not _check_matching_var_meta(datasets):
        message = (
            "Variable meta attributes vary between different tajectories. "
            "This indicates inconsistencies like distinct units between trajectories. "
            "Please ensure consistency between datasets before merging."
        )
        logging.warning(f"{message} Merge result may be inconsistent.")
        # TODO: Do we want to merge anyway?
        raise ValueError(f"{message} Will not merge.")

    consistent_metadata, distinct_metadata = _merge_traj_metadata(datasets)

    trajids = distinct_metadata["trajid"]

    datasets = [ds.expand_dims(trajid=[id]) for ds, id in zip(datasets, trajids)]

    # trajids = pd.Index(meta["trajid"], name="trajid")
    # coords_trajids = xr.Coordinates(indexes={"trajid": trajids})
    # breakpoint()
    layers = xr.concat(datasets, dim="trajid", combine_attrs="drop_conflicts")

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
            k: ("trajid", v, {"description:": f"Attribute {k} merged in concatenation"})
            for k, v in distinct_metadata
            if k != "trajid"
        }
    )

    layers.attrs["is_multi_trajectory"] = True
    if not isinstance(layers, Trajectory):
        layers = Trajectory(layers)

    if TYPE_CHECKING:
        assert isinstance(layers, Trajectory)

    return layers
