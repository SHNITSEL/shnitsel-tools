from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm

from shnitsel.data.tree.child_support_functions import find_child_key

from . import (
    ShnitselDBRoot,
    CompoundGroup,
    DataGroup,
    DataLeaf,
    TreeNode,
)


DataType = TypeVar("DataType")


def build_shnitsel_db(
    data: ShnitselDBRoot[DataType]
    | CompoundGroup[DataType]
    | DataGroup[DataType]
    | DataLeaf[DataType]
    | TreeNode[Any, DataType]
    | DataType
    | Sequence[CompoundGroup[DataType]]
    | Sequence[DataGroup[DataType] | DataLeaf[DataType]]
    | Sequence[TreeNode[Any, DataType]]
    | Sequence[DataType],
    dtype: type[DataType] | UnionType | None = None,
) -> ShnitselDBRoot[DataType]:
    """
    Function to generate a full -- i.e. up to ShnitselDBRoot -- Shnitsel DB structure.

    Wraps trajectories in DataLeaf structures, extends the tree with missing parent structures.

    Parameters
    ----------
    data : ShnitselDBRoot[DataType] | CompoundGroup[DataType] | DataGroup[DataType] | DataLeaf[DataType] | DataType | Sequence[CompoundGroup[DataType]] | Sequence[DataGroup[DataType]  |  DataLeaf[DataType]] | Sequence[DataType]
        Input data to be wrapped in a ShnitselDB format
    dtype :  type[DataType] | UnionType, optional
        The datatype that data in this tree should have.
        If not provided will be inferred from the data in the tree.

    Returns
    -------
    ShnitselDBRoot[DataType]
        The resulting ShnitselDB dataset structure.

    Raises
    ------
        ValueError: If an unsupported `data` argument was provided.
        ValueError: If a list of xr.DataTree objects on incompatible Levels of the ShnitselDB hierarchy was provided, e.g. a mix of DataLeaf and CompoundGroup nodes.
        ValueError: If the provided data is of no ShnitselDB format type.
    """
    if isinstance(data, Sequence):
        is_likely_root = all(isinstance(child, CompoundGroup) for child in data)
        if is_likely_root:
            try:
                return ShnitselDBRoot(
                    compounds={str(child.name): child for child in data},  # type: ignore # The above check ensures the children to be Compound group
                    dtype=dtype,
                )
            except:
                pass

        is_likely_compound = all(
            (isinstance(child, DataGroup) and not isinstance(child, CompoundGroup))
            or isinstance(child, DataLeaf)
            for child in data
        )
        if is_likely_compound:
            return build_shnitsel_db(
                CompoundGroup(
                    children={str(child.name): child for child in data},  # type: ignore # The above check ensures the children to be DataGroup or DataLeaf instances
                    dtype=dtype,
                ),
                dtype=dtype,
            )

        # We have raw data in here, wrap it:
        res: dict[str, DataGroup[DataType] | DataLeaf[DataType]] = {}
        for child in data:
            if isinstance(child, (DataGroup, DataLeaf)):
                res[find_child_key(res.keys(), child, "data")]
            elif isinstance(child, TreeNode):
                raise ValueError(
                    "Unsupported child type at data level: %s" % type(child)
                )
            else:
                name_candidate: str | None = None
                if name_candidate is None:
                    name_candidate = getattr(child, 'name', None)
                if name_candidate is None:
                    name_candidate = getattr(child, 'trajid', None)
                if name_candidate is None:
                    name_candidate = getattr(child, 'trajectory_id', None)
                if name_candidate is None:
                    name_candidate = getattr(child, 'id', None)

                new_child = DataLeaf[DataType](
                    name=name_candidate, data=child, dtype=dtype
                )
                child_key = find_child_key(res.keys(), new_child, "data")
                new_child._name = child_key
                res[child_key] = new_child
        return build_shnitsel_db(
            CompoundGroup(
                children={str(key): child for key, child in res.items()}, dtype=dtype
            ),
            dtype=dtype,
        )  # type: ignore # The above check ensures the children to be Compound groups
    if isinstance(data, ShnitselDBRoot):
        # If we already are a root node, return the structure.
        return data
    elif isinstance(data, CompoundGroup):
        key = find_child_key([], data, "compound")
        return ShnitselDBRoot(compounds={key: data}, dtype=dtype)
    elif isinstance(data, DataGroup):
        key = find_child_key([], data, "group")
        return build_shnitsel_db(
            CompoundGroup(children={key: data}, dtype=dtype), dtype=dtype
        )
    elif isinstance(data, DataLeaf):
        key = find_child_key([], data, "data")
        return build_shnitsel_db(
            CompoundGroup(children={key: data}, dtype=dtype), dtype=dtype
        )
    else:
        raise ValueError("Unsupported node type provided: %s" % type(data))


complete_shnitsel_tree = build_shnitsel_db
