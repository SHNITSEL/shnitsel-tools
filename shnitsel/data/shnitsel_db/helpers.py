import logging
from typing import Iterable, Mapping

from .db_trajectory_data import TrajectoryData
from .db_trajectory_group import TrajectoryGroup


def traj_list_to_child_mapping(
    list: Iterable[TrajectoryData | TrajectoryGroup],
) -> Mapping[str, TrajectoryData | TrajectoryGroup]:
    """Function to set the labels in a children-dict for a DataTree.

    Will enumerate individual trajectories and assign the group name to groups.

    Args:
        list (Iterable[TrajectoryData  |  TrajectoryGroup]): The list of groups and Trajectories

    Returns:
        Mapping[str, TrajectoryData | TrajectoryGroup]: The dict mapping labels to children
    """
    # TODO: FIXME: Deal with conflicting group names
    res = {}

    for i, v in enumerate(list):
        if isinstance(v, TrajectoryData):
            if v.dataset is not None and "trajid" in v.dataset.attrs:
                trajid = str(v.dataset.attrs["trajid"])
            else:
                trajid = f"{i}"

            if trajid in res:
                trajid = "_" + trajid

            res = {trajid: v}
        elif isinstance(v, TrajectoryGroup):
            key = v.name if v.name is not None else f"{i}"

            if key in res:
                key = "_" + key

            res = {key: v}
        else:
            logging.error(
                f"Invalid type: {type(v)}, only TrajectoryData or TrajectoryGroup allowed."
            )
    return res
