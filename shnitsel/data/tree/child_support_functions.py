from typing import Hashable, Iterable
from .node import TreeNode


def find_child_key(
    existing_child_keys: Iterable[Hashable],
    new_child: TreeNode,
    default_prefix: str = "child",
) -> str:
    """Helper function to find a new collision free key for a new child node
    given the existing set of keys and the new node

    Will first attempt to extract a name candidate from the node and then
    try to resolve potential name collisions by appending suffixes.

    Parameters
    ----------
    existing_child_keys : Iterable[Hashable]
        Set of existing keys in use for existing children.
    new_child : TreeNode
        The new node for which we want to find a key.
    default_prefix : str, optional
        A default prefix to use if we cannot extract a name candidate, by default "child"

    Returns
    -------
    str
        A derived key that does not have any collisions in `existing_child_keys`
        and incorporates information from `new_child` if possible.

    Raises
    ------
    OverflowError
        If no unused child name has been found after 1000 attempts.
    """
    from .data_leaf import DataLeaf

    try:
        name_candidate = new_child.name
    except:
        name_candidate = None

    if name_candidate is None:
        if isinstance(new_child, DataLeaf) and new_child.has_data:
            name_candidate = get_data_name_candidate(new_child.data)

    if name_candidate is None:
        name_candidate = default_prefix
    elif name_candidate not in existing_child_keys:
        return name_candidate

    for i in range(1000):
        new_candidate = name_candidate + f"_{i}"

        if new_candidate not in existing_child_keys:
            return new_candidate

    raise OverflowError(
        "Could not find a new child name that was not already in use in `1000` attempts."
    )


def get_data_name_candidate(data, name_candidate: str | None = None) -> str | None:
    """Helper function to extract name candidates from arbitrary data
    that could be stored in a Tree leave.

    Parameters
    ----------
    data : Any
        The arbitrary data for which we want to derive a name
    name_candidate : str | None, optional
        A potentially already existing name candidate that may be used, by default None

    Returns
    -------
    str | None
        A name based on the `data` or `None` if no name could be derived.
    """
    if name_candidate is None:
        name_candidate = getattr(data, 'name', None)
    if name_candidate is None:
        name_candidate = getattr(data, 'trajid', None)
    if name_candidate is None:
        name_candidate = getattr(data, 'id', None)
    return name_candidate
