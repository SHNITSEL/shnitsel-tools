from typing import Any, Sequence, TypeVar

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
) -> ShnitselDBRoot[DataType]:
    """
    Function to generate a full -- i.e. up to ShnitselDBRoot -- Shnitsel DB structure.

    Wraps trajectories in DataLeaf structures, extends the tree with missing parent structures.

    Parameters
    ----------
    data : ShnitselDBRoot[DataType] | CompoundGroup[DataType] | DataGroup[DataType] | DataLeaf[DataType] | DataType | Sequence[CompoundGroup[DataType]] | Sequence[DataGroup[DataType]  |  DataLeaf[DataType]] | Sequence[DataType]
        Input data to be wrapped in a ShnitselDB format

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
                    compounds={str(child.name): child for child in data}  # type: ignore # The above check ensures the children to be Compound groups
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
                CompoundGroup(children={str(child.name): child for child in data})  # type: ignore # The above check ensures the children to be DataGroup or DataLeaf instances
            )

        # We have raw data in here, wrap it:
        res: Sequence[DataGroup[DataType] | DataLeaf[DataType]] = []
        for child in data:
            if isinstance(child, (DataGroup, DataLeaf)):
                res.append(child)
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
                    name_candidate = getattr(child, 'id', None)

                res.append(DataLeaf[DataType](name=name_candidate, data=child))
            return build_shnitsel_db(
                CompoundGroup(children={str(child.name): child for child in res})
            )  # type: ignore # The above check ensures the children to be Compound groups
    if isinstance(data, ShnitselDBRoot):
        # If we already are a root node, return the structure.
        return data
    elif isinstance(data, CompoundGroup):
        # If we are provided a single trajectory wrap in a TrajectoryData node and build rest of tree
        return ShnitselDBRoot(compounds={data.name: data})
    elif isinstance(data, DataGroup):
        # If we are provided a single trajectory wrap in a TrajectoryData node and build rest of tree
        return build_shnitsel_db(CompoundGroup(children={data.name: data}))
    elif isinstance(data, DataLeaf):
        return build_shnitsel_db(CompoundGroup(children={data.name: data}))
    else:
        raise ValueError("Unsupported node type provided: %s" % type(data))


complete_shnitsel_tree = build_shnitsel_db
