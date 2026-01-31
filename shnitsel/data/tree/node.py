import abc
from collections.abc import Iterable
from dataclasses import dataclass
from functools import update_wrapper, wraps
import logging
import os
from types import GenericAlias, UnionType
from typing import (
    Any,
    Callable,
    Hashable,
    Iterator,
    Mapping,
    Self,
    Sequence,
    TypeVar,
    Generic,
    overload,
)
import typing
from typing_extensions import Literal, TypeForm

from shnitsel.core.typedefs import ErrorOptionsWithWarn
from shnitsel.data.dataset_containers.multi_layered import MultiSeriesLayered
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from ..trajectory_grouping_params import TrajectoryGroupingMetadata
from pathlib import Path


ChildType = TypeVar("ChildType", bound="TreeNode|None", covariant=True)
DataType = TypeVar("DataType", covariant=True)
NewDataType = TypeVar("NewDataType")
NewChildType = TypeVar("NewChildType", bound="TreeNode|None")
ResType = TypeVar("ResType", covariant=True)
KeyType = TypeVar("KeyType")

T = TypeVar("T")

_class_cache = {}


@dataclass
class TreeNode(Generic[ChildType, DataType], abc.ABC):
    """Base class to model a tree structure of arbitrary data type to keep
    trajectory data with hierarchical structure in.

    Has two type parameters to allow for explicit type checks:
    - `ChildType`: Which node types are allowed to be registered as children of this node.
    - `DataType`: What kind of data is expected within this tree if the data is not None.
    """

    @classmethod
    def _get_extended_class_name(cls: type, datatypes: Sequence[type]) -> str:
        dtype_string = "|".join([ot.__name__ for ot in datatypes])
        resname = f"{cls.__name__}[{dtype_string}]"
        return resname

    @classmethod
    def _create_extended_node_class(
        cls: type[Self], datatypes: list[tuple[type, list[str], list[str]]]
    ) -> type[Self]:
        """Create a new version of the class with added methods for the datatypes."""

        def make_mapped_method(method_name: str, docstring: str | None = None):
            # def method(self, *args, **kw):
            #     return getattr(object.__getattribute__(self, "_obj"), name)(*args, **kw)
            def method_wrapper(self: TreeNode, *args, **kwargs):
                # TODO: FIXME: deal with trees in the args?
                def internally_mapped_method(data, *iargs, **ikwargs):
                    # There should never be arguments set directly on this method call.
                    if hasattr(data, method_name):
                        tmp_res = getattr(data, method_name)
                        if callable(tmp_res):
                            # We do not want callable types in the tree, that would be weird
                            return tmp_res(*args, **kwargs)
                        else:
                            logging.warning(
                                "Access of method %s on data at `%s` is clashing with non-callable property. The result will be skipped.",
                                method_name,
                                self.path,
                            )
                    return None

                return self.map_data(internally_mapped_method)
                # return getattr(object.__getattribute__(self, "_obj"), name)(*args, **kw)

            method_wrapper.__doc__ = docstring
            return method_wrapper

        def make_mapped_property(prop_name: str, docstring: str | None = None):
            def prop_wrapper(self: TreeNode):
                def acessor_wrapper(data):
                    if hasattr(data, prop_name):
                        tmp_res = getattr(data, prop_name)
                        if not callable(tmp_res):
                            # We do not want callable types in the tree, that would be weird
                            return tmp_res
                        else:
                            logging.error(
                                "Access of property %s would yield callable entries in the tree at `%s`, which is not supported. The result will be skipped.",
                                prop_name,
                                self.path,
                            )
                    return None

                return self.map_data(acessor_wrapper)
                # return getattr(object.__getattribute__(self, "_obj"), name)(*args, **kw)

            return property(fget=prop_wrapper, doc=docstring)

        namespace = {}
        for key in dir(cls):
            # Copy base class properties
            namespace[key] = getattr(cls, key, None)
            if namespace[key] is None:
                del namespace[key]

        # Patch in data class wrapper accessors
        for otype, props, methods in datatypes:
            # if hasattr(theclass, name):
            # Do not override things we explicitly override

            for prop in props:
                if prop in namespace:
                    # Do not override existing entries in namespace
                    continue
                prop_entry = getattr(otype, prop, None)
                docstring = None
                if prop_entry is not None:
                    if hasattr(prop_entry, "__doc__"):
                        docstring = getattr(prop_entry, "__doc__")
                namespace[prop] = make_mapped_property(prop, docstring=docstring)

            for method in methods:
                if method in namespace:
                    # Do not override existing entries in namespace
                    continue
                method_entry = getattr(otype, method, None)
                if callable(method_entry):
                    docstring = None
                    if method_entry is not None:
                        if hasattr(method_entry, "__doc__"):
                            docstring = getattr(method_entry, "__doc__")
                    namespace[method] = update_wrapper(
                        make_mapped_method(method, docstring=docstring), method_entry
                    )

            # if hasattr(theclass, name) and not hasattr(cls, name):
            #     namespace[name] = make_method(name)
        configured_datatypes = [ot for ot, _, _ in datatypes]
        namespace["__configured_datatypes__"] = configured_datatypes
        resname = cls._get_extended_class_name(configured_datatypes)
        logging.debug(f"Creating patched node class {resname}")

        # return type("%s(%s)" % (cls.__name__, theclass.__name__), (cls,), namespace)
        return type(resname, (cls,), namespace)

    def __class_getitem__(
        cls: type[Self], args: "TypeVar | tuple[TypeVar , ...]"
    ) -> type[Self]:
        # print(typing.get_args(cls))
        # print(typing.get_args(args))
        # print(f"{cls=}[{args=}]")
        data_type_options = None
        if not isinstance(args, TypeVar):
            if isinstance(args, tuple):
                datatype_arg = args[-1]
            else:
                datatype_arg = args
        else:
            datatype_arg = None

        base = super().__class_getitem__(args)
        # print(base, cls)
        # if not issubclass(base, cls):
        #     base = cls

        if datatype_arg is not None:
            if not isinstance(datatype_arg, TypeVar):
                if isinstance(datatype_arg, type):
                    data_type_options = [datatype_arg]
                elif isinstance(datatype_arg, UnionType):
                    data_typeSelf_options = list(typing.get_args(datatype_arg))

                # print(f"{part_types=}")

        if data_type_options is not None:
            resname = cls._get_extended_class_name([ot for ot in data_type_options])
            class_cache = getattr(cls, "__class_cache__", {})
            if not resname in class_cache:
                types_with_props_and_methods = []
                for otype in data_type_options:
                    # Find non-private members
                    non_private_entries = [
                        d for d in dir(otype) if not d.startswith("_")
                    ]

                    public_properties = [
                        d
                        for d in non_private_entries
                        if not callable(getattr(otype, d, None))
                    ]
                    public_funcs = [
                        d
                        for d in non_private_entries
                        if callable(getattr(otype, d, None))
                    ]

                    types_with_props_and_methods.append(
                        (otype, public_properties, public_funcs)
                    )
                logging.debug(
                    "Creating new class for dtype fields and methods: %s",
                    types_with_props_and_methods,
                )
                resclass = cls._create_extended_node_class(types_with_props_and_methods)
                class_cache[resname] = resclass
                setattr(cls, "__class_cache__", class_cache)
            return class_cache[resname]
        return base

    _name: str | None
    _dtype: type[DataType] | UnionType | None

    _data: DataType | None
    _children: Mapping[Hashable, ChildType]
    _attrs: Mapping[str, Any]

    _parent: Self | None
    _level_name: str | None

    def __new__(
        cls: type[Self],
        *,
        data: DataType | None = None,
        children: Mapping[Hashable, ChildType] | None = None,
        dtype: type[DataType] | UnionType | None = None,
        **kwargs,
    ) -> Self:
        if dtype is not None:
            ndtype = typing.get_origin(dtype)
            if ndtype is not None:
                dtype = ndtype
            assert dtype is not None

            filled_in_dtype = dtype
            if data is not None:
                if isinstance(dtype, UnionType):
                    assert isinstance(data, typing.get_args(dtype)), (
                        "Provided data did not match provided dtype"
                    )
                else:
                    assert isinstance(data, dtype), (
                        "Provided data did not match provided dtype"
                    )
        else:
            if data is not None:
                # If we have data, try and use the type of that data
                filled_in_dtype = type(data) or None
            else:
                filled_in_dtype = TreeNode._dtype_guess_from_children(children)

        dtype = filled_in_dtype
        if data is not None:
            kwargs['data'] = data
        if children is not None:
            kwargs['children'] = children
        if dtype is not None:
            kwargs['dtype'] = dtype

        if isinstance(cls, GenericAlias):
            cls = typing.get_origin(cls)

        base_class: type[Self] = cls

        if not hasattr(cls, "__configured_datatypes__") and dtype is not None:
            if isinstance(cls, GenericAlias):
                orig_class = typing.get_origin(cls)
                cls = orig_class or cls

            base_class = cls[dtype]

        logging.debug("Creating object for %s", base_class)
        try:
            obj = object.__new__(base_class)
        except TypeError:
            logging.error(
                "Failed to create object of type: %s because of generic type settings. You may have specified `dtype` parameters with generic aliases. (class parameter: %s)",
                base_class,
                cls,
            )
            try:
                obj = super().__new__(base_class)
            except TypeError as e:
                # print(e)
                # print(base_class)
                logging.error(
                    "Failed to create object of type: %s because of generic type settings. You may have specified `dtype` parameters with generic aliases. (class parameter: %s)",
                    base_class,
                    cls,
                )
                raise TypeError(
                    "Failed to create object of type: %s because of generic type settings. You may have specified `dtype` parameters with generic aliases. (class parameter: %s)"
                    % (base_class, cls)
                ) from e

        obj.__init__(**kwargs)

        return obj

    def __init__(
        self,
        *,
        name: str | None,
        data: DataType | None = None,
        children: Mapping[Hashable, ChildType] | None = None,
        attrs: Mapping[str, Any] | None = None,
        level_name: str | None = None,
        dtype: type[DataType] | None = None,
        **kwargs,
    ):
        self._name = name
        self._data = data
        self._children = children if children is not None else dict()
        self._attrs = attrs if attrs is not None else dict()
        self._parent = None
        self._level_name = (
            level_name if level_name is not None else self.__class__.__qualname__
        )
        # Set parent entry
        for child in self._children.values():
            if child is not None:
                child._parent = self

        if dtype is not None:
            ndtype = typing.get_origin(dtype)
            if ndtype is not None:
                dtype = ndtype
            assert dtype is not None

            filled_in_dtype = dtype
            if data is not None:
                if isinstance(dtype, UnionType):
                    assert isinstance(data, typing.get_args(dtype)), (
                        "Provided data did not match provided dtype"
                    )
                else:
                    assert isinstance(data, dtype), (
                        "Provided data did not match provided dtype"
                    )
        else:
            if data is not None:
                # If we have data, try and use the type of that data
                filled_in_dtype = type(data)
            else:
                filled_in_dtype = TreeNode._dtype_guess_from_children(self._children)

        self._dtype = filled_in_dtype

    @staticmethod
    def _dtype_guess_from_children(
        children: Mapping | None,
    ) -> type | UnionType | None:
        if children:
            # If we have children, try and construct the type from the basis type
            child_types_collection: set[type | UnionType] = set()

            for child in children.values():
                if child is not None:
                    cdtype = child.dtype
                    if cdtype is not None:
                        child_types_collection.add(cdtype)

            filled_in_dtype = None
            for ct in child_types_collection:
                if filled_in_dtype is None:
                    filled_in_dtype = ct
                else:
                    filled_in_dtype = filled_in_dtype | ct
                    # print(filled_in_dtype)
        else:
            # If nothing helps, don't add a dtype for this node
            filled_in_dtype = None

        return filled_in_dtype

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, ChildType] | None = None,
        dtype: None = None,
        data: DataType | None = None,
        **kwargs,
    ) -> Self: ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: None = None,
        **kwargs,
    ) -> "TreeNode[NewChildType, ResType]": ...

    @overload
    def construct_copy(
        self,
        children: None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> "TreeNode[Any, ResType]": ...

    @abc.abstractmethod
    def construct_copy(
        self,
        children: Mapping[Hashable, ChildType]
        | Mapping[Hashable, NewChildType]
        | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> Self | "TreeNode[NewChildType, ResType]":
        """Every class inheriting from TreeNode should implement this method to create a copy of that subtree
        with appropriate typing or just plain up creating a copy of the subtree, if no updates are requested.

        Support for changing the typing by changing child types, setting the explicit `dtype` or by providing
        a new `data` entry should be supported by the base class.

        Parameters
        ----------
        data : ResType | None, optional
            The new data to be set in the copy of this node, by default None, which should populate it with the node's current data
        children : Mapping[str, NewChildType], optional
            A new set of children to replace the old mapping of children can be provided with this parameter.
            The data type can also be changed with appropriate typing here:
        dtype : type[ResType] | UnionType | None, optional
            An explicit argument to set the `dtype` property of the new subtree, by default None.

        Returns
        -------
        Self | TreeNode[TreeNode[Any, RestType]|None, ResType]
            Returns a new subtree with a duplicate of this node in regards to metadata at its root and
            updates properties as provided.
        """
        raise NotImplementedError(
            "The node type should provide a specialization of the `construct_copy()` function"
        )

    # @overload
    # def construct_copy(
    #     self,
    #     data: ResType,
    #     dtype: None = None,
    #     **kwargs,
    # ) -> "TreeNode[TreeNode[Any, ResType]|None, ResType]": ...

    # @overload
    # def construct_copy(
    #     self,
    #     data: None,
    #     dtype: type[ResType] | TypeForm[ResType],
    #     **kwargs,
    # ) -> "TreeNode[TreeNode[Any, ResType]|None, ResType]": ...

    # @overload
    # def construct_copy(
    #     self,
    #     data: DataType | None,
    #     dtype: type[DataType] | TypeForm[DataType] | None,
    #     **kwargs,
    # ) -> Self: ...

    # def construct_copy(
    #     self,
    #     data: ResType | None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     **kwargs,
    # ) -> Self | "TreeNode[Any, ResType]":
    #     """Function to construct a copy with optionally a new `dtype`
    #     set for this node and the data value updated.

    #     Parameters
    #     ----------
    #     data: DataType | ResType, optional
    #         New data for the copy node. Either with the same datatype as the current tree or
    #         with a new ResType type such that a new tree with a different `DataType` template argument
    #         can be constructed
    #     dtype: type[DataType] | TypeForm[DataType] | type[ResType] | TypeForm[ResType]
    #         The data type of the data in this tree.
    #     **kwargs: dict[str, Any], optional
    #         A dictionary holding the settings for the constructor of the current node's type.

    #     Returns
    #     -------
    #     Self
    #         Either the same type as self if data and dtype have not been updated
    #     TreeNode[TreeNode[Any, ResType], ResType]
    #         Or effectively the same original type with update DataType variable
    #     """
    #     if data is not None:
    #         kwargs['data'] = data
    #     elif 'data' not in kwargs:
    #         kwargs['data'] = self._data

    #     if dtype is not None:
    #         kwargs['dtype'] = dtype
    #     elif 'dtype' not in kwargs:
    #         kwargs['dtype'] = self._dtype

    #     if 'name' not in kwargs:
    #         kwargs['name'] = self._name
    #     if 'children' not in kwargs:
    #         kwargs['children'] = {
    #             # TODO: FIXME: Figure out this typing issue
    #             k: v.construct_copy(dtype=dtype)
    #             for k, v in self._children.items()
    #             if v is not None
    #         }
    #     if 'attrs' not in kwargs:
    #         kwargs['attrs'] = dict(self._attrs)
    #     if 'level_name' not in kwargs:
    #         kwargs['level_name'] = str(self._level_name)

    #     return type(self)(
    #         **kwargs,
    #     )

    @property
    def path(self) -> str:
        if self._parent is None:
            return "/"
        else:
            parent_key = [k for k, v in self._parent._children.items() if v == self]
            if parent_key:
                parent_key = parent_key[0]
            else:
                "?"
            return self._parent.path + str(parent_key) + "/"

    def __len__(self) -> int:
        """Returns the `size` of this node, i.e. how many children it has.

        Be aware that this means that it will return 0 for Leaf nodes that may hold data.

        Returns
        -------
        int
            The number of children of this node
        """
        return len(self._children)

    def __contains__(self, value: str | ChildType) -> bool:
        if isinstance(value, str):
            return value in self._children
        else:
            return value in self._children.values()

    def __getitem__(
        self, key: str | tuple[str]
    ) -> "TreeNode[Any, DataType] | DataType | None":
        path_parts: tuple
        if isinstance(key, str):
            path_parts = Path(key).parts
        elif isinstance(key, tuple):
            path_parts = key
        else:
            raise ValueError("Unsupported index type: %s", type(key))

        if len(path_parts) == 0:
            return self

        first_part = path_parts[0]
        path_tail = tuple(path_parts[1:])
        if first_part == os.path.sep:
            # print("goto root")
            return self.root[path_tail]
        elif first_part == '.':
            # print("goto self")
            return self[path_tail]
        elif first_part == '..':
            if self._parent is not None:
                # print("goto parent")
                return self._parent[path_tail]
            # print("goto parent impossible")
            return self[path_tail]
        else:
            if self.has_data and first_part == 'data':
                # print("yield data")
                return self.data
            elif first_part in self._children:
                child_entry = self._children[first_part]
                if child_entry is None:
                    return None
                else:
                    return child_entry[path_tail]
        return None

    def __setitem__(self, key, value):
        # self.daten[key] = value
        raise AssertionError(
            "Cannot set tree items with the array syntax. Would violate the "
        )

    @property
    def is_leaf(self) -> bool:
        return len(self._children) == 0

    @property
    def has_data(self) -> bool:
        return self._data is not None

    @property
    def dtype(self) -> type[DataType] | UnionType | None:
        return self._dtype

    @property
    def data(self) -> DataType:
        if self._data is None:
            raise ValueError("Object has no data to be retrieved")
        else:
            return self._data

    @property
    def children(self) -> Mapping[Hashable, ChildType]:
        return self._children

    @property
    def root(self) -> "TreeNode[Any, DataType]":
        if self._parent is None:
            return self
        else:
            return self._parent.root

    @property
    def attrs(self) -> Mapping[str, Any]:
        return self._attrs

    @property
    def name(self) -> str:
        if self._name is None:
            raise KeyError("Node has no `name` attribute set.")
        else:
            return self._name

    # def map_over_child_nodes(
    #     self, func: Callable[[Self], ResType | None]
    # ) -> Mapping[Hashable, ResType]:
    #     new_children = {
    #         k: res
    #         for k, v in self._children.items()
    #         if v is not None and (res := v.map_node(func)) is not None
    #     }
    #     return new_children

    def map_subtree(self, func: Callable[[Self], ResType]) -> ResType:
        """Just a helper function with telling name to apply a function
        to the root node of this current subtree.

        Simply calls `func(self)`.

        Parameters
        ----------
        func : Callable[[Self], ResType]
            The function to apply to this node

        Returns
        -------
        ResType
            The result of `funct(self)`.
        """
        return func(self)

    @abc.abstractmethod
    def group_children_by(
        self,
        key_func: Callable[["TreeNode"], KeyType],
        group_leaves_only: bool = False,
    ) -> Self | None:
        """Method to group nodes within this current subtree by keys
        as retrieved via `key_func`.

        Can be used to group data within this tree by metadata, e.g.
        to separate trajectory data with different simulation settings into
        distinct groups.

        Adds new groups into the tree structure.

        Parameters
        ----------
        key_func : Callable[[TreeNode], KeyType]
            Key function that should map Any tree node that is not excluded, e.g. by setting
            `group_leaves_only` to a key value that should be a dataclass and should be
            equal for two nodes if and only if those nodes should eventually end up in the same group.
        group_leaves_only : bool, optional
            Flag to control whether grouping should only be applied to
            `DataLeaf` nodes, by default False

        Returns
        -------
        Self | None
            The current node after its subtree has been grouped.
            If no keys could be retrieved, the result may be `None`.
        """
        ...

    @overload
    def map_data(
        self,
        func: Callable[..., ResType | None] | Callable[..., ResType],
        *args,
        keep_empty_branches: Literal[True] = True,
        dtype: type[ResType],
        **kwargs,
    ) -> "TreeNode[Any,ResType]": ...

    @overload
    def map_data(
        self,
        func: Callable[..., ResType | None] | Callable[..., ResType],
        *args,
        keep_empty_branches: Literal[False],
        dtype: type[ResType],
        **kwargs,
    ) -> "TreeNode[Any,ResType]|None": ...

    @overload
    def map_data(
        self,
        func: Callable[..., ResType | None] | Callable[..., ResType],
        *args,
        keep_empty_branches: Literal[False],
        dtype: None = None,
        **kwargs,
    ) -> "TreeNode[Any,ResType]|None": ...

    @overload
    def map_data(
        self,
        func: Callable[..., ResType | None] | Callable[..., ResType],
        *args,
        keep_empty_branches: Literal[True] = True,
        dtype: None = None,
        **kwargs,
    ) -> "TreeNode[Any,ResType]": ...

    @overload
    def map_data(
        self,
        func: Callable,
        *args,
        keep_empty_branches: Literal[False],
        dtype: None = None,
        **kwargs,
    ) -> "TreeNode|None": ...

    @overload
    def map_data(
        self,
        func: Callable,
        *args,
        keep_empty_branches: Literal[True] = True,
        dtype: None = None,
        **kwargs,
    ) -> "TreeNode": ...

    def map_data(
        self,
        func: Callable[..., ResType | None],
        *args,
        keep_empty_branches: bool = True,
        dtype: type[ResType] | None = None,
        **kwargs,
    ) -> "TreeNode[Any,ResType]|TreeNode|None":
        """Helper function to apply a mapping function to all data in leaves of this tree

        The function `func` is applied to all `DataLeaf` instances with `data` within them.
        If `keep_empty_branches=False` is set, will truncate branches without any data in them or without any further children.

        Parameters
        ----------
        func : Callable[[DataType], ResType  |  None]
            The mapping function to apply to data in this subtree.
        keep_empty_branches : bool, optional
            Flag to control whether branches/subtrees without any data in them should be truncated, by default False to keep the same structure
        dtype : type[ResType] | None, optional
            Optional parameter to explicitly specify the `dtype` for the resulting tree, by default None
        *args
            Positional arguments to pass to the call to `func`
        **kwargs
            Keyword-style arguments to pass to the call to `func`

        Returns
        -------
        TreeNode[Any,ResType]|None
            The resulting node after the subtree has been mapped or None if truncation is active and the subtree has no data after mapping.
        """
        if not self.is_leaf:
            new_children: dict[Hashable, TreeNode[Any, ResType]] = {
                k: res
                for k, v in self._children.items()
                if v is not None
                and (
                    res := v.map_data(
                        func,
                        *args,
                        keep_empty_branches=keep_empty_branches,
                        dtype=dtype,
                        **kwargs,
                    )
                )
                is not None
            }

            if len(new_children) == 0 and not keep_empty_branches:
                return None
            else:
                return self.construct_copy(children=new_children, dtype=dtype)
        else:
            # Map data in leaves
            if self.has_data:
                new_data = func(self.data, *args, **kwargs)
            else:
                new_data = None

            if not keep_empty_branches and new_data is None:
                return None
            else:
                # This yields a different kind of tree.
                return self.construct_copy(data=new_data, dtype=dtype)

    def map_filtered_nodes(
        self,
        filter_func: Callable[["TreeNode[Any, DataType]"], bool],
        map_func: Callable[["TreeNode[Any, DataType]"], "TreeNode[Any, ResType]|None"],
        dtype: type[ResType] | None = None,
    ) -> "TreeNode[Any, ResType]|None":
        """Map nodes using `map_func()` if the filter function `filter_func` picks them as relevant.

        If the node is not picked by `filter_func` a copy will be created with its children being recursively mapped
        according to the same rule.
        If a node is mapped, the mapping function `map_func` must take care of potential mapping over children.

        Parameters
        ----------
        filter_func : Callable[[TreeNode[Any, DataType]], bool]
            Filter function to apply to nodes in the current subtree of any kind. Must return `True` for all nodes to which `map_func` should be applied.
        map_func : Callable[[TreeNode[Any, DataType]], TreeNode[Any, ResType]|None]
            Mapping function that transforms a selected node of a certain datatype to a consistent new data type `RestType`.
        dtype : type[ResType] | None, optional
            Optional parameter to explicitly specify the `dtype` for the resulting tree, by default None.

        Returns
        -------
        TreeNode[Any, ResType]
            A new subtree with the data type changed and select subtrees mapped.
        None
            If the node was filtered and the map function returned None
        """
        # from .data_group import DataGroup
        from .data_leaf import DataLeaf
        # from .compound import CompoundGroup
        # from .tree import ShnitselDBRoot

        new_node: TreeNode[Any, ResType] | None

        if filter_func(self):
            new_node = map_func(self)
        else:
            if isinstance(self, DataLeaf):
                new_node = self.construct_copy(
                    data=self._data,  # type: ignore # If the user gives us a conflicting type that does not include unmapped leaves, it is their problem.
                    dtype=dtype,  # type: ignore # Same as in the prior line
                )
            else:
                # assert isinstance(self, (DataGroup, CompoundGroup, ShnitselDBRoot)), (
                #     "Unsupported node type provided to `map` function: %s" % type(self)
                # )

                new_children = {
                    k: res
                    for k, v in self.children.items()
                    if v is not None
                    and (
                        res := v.map_filtered_nodes(filter_func, map_func, dtype=dtype)
                    )
                    is not None
                }
                new_node = self.construct_copy(
                    children=new_children,
                    dtype=dtype,  # type: ignore # By mapping the children, we are sure that they now also hold the resulting datatype.
                )

        return new_node

    def filter_nodes(
        self,
        filter_func: Callable[..., bool],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> Self | None:
        """Function to filter the nodes in this tree and create a new tree that are ancestors of
        at least one accepted node.

        If `keep_empty_branches=True`, all branches in which there are no accepted nodes, will be truncated.
        If `filter_func` does not return `True`, the entire subtree starting at this node, will be dropped.

        Parameters
        ----------
        filter_func : Callable[..., bool]
            A filter function that should return True for Nodes that should be kept within the Tree and `False` for Nodes that should be kicked out together with their entire subtree.
        recurse : bool, optional
            Whether to recurse the filtering into the children of kept nodes, by default True
        keep_empty_branches : bool, optional
            A flag to enable truncation of branches with only empty lists of children and no data, by default False

        Returns
        -------
        Self | None
            Either a copy of the current subtree if it is kept or None if the subtree is omitted
        """
        from .data_group import DataGroup
        from .data_leaf import DataLeaf
        from .compound import CompoundGroup
        from .tree import ShnitselDBRoot

        assert isinstance(self, (DataGroup, CompoundGroup, ShnitselDBRoot)), (
            "Unsupported node type provided to `map` function: %s" % type(self)
        )

        keep_self = filter_func(self)

        # Stop if the node is not kept.
        if not keep_self:
            return None

        if isinstance(self, DataLeaf):
            return self.construct_copy()
        else:
            new_children = None
            if recurse:
                new_children = {
                    k: res
                    for k, v in self._children.items()
                    if v is not None
                    and (res := v.filter_nodes(filter_func)) is not None
                }

                if len(new_children) == 0:
                    new_children = None

            if not keep_empty_branches and not keep_self and new_children is None:
                return None
            else:
                tmp_res: Self = self.construct_copy(children=new_children)  # type: ignore # Our filtering effectively only copies the subtree, no type modification is performed.
                return tmp_res

    def add_child(self, child_name: str | None, child: ChildType) -> Self:
        """Add a new child node with a preferred name in the mapping of children.
        If the child name is already in use, will attempt to find a collision-free alternative name.

        Parameters
        -----------
            child_name (str | None): The preferred name under which the child should be registered.
                    To avoid overriding, a different name will be chosen if the name is in use.
            child (ChildType): Object to register as the child-subtree

        Raises
        -----------
            OverflowError: If the attempts to find a new collision-free name have exceeded 1000.

        Returns
        -----------
            Self: The new instance of a subtree
        """
        new_children = dict(self._children)
        if child_name is not None and child_name not in new_children:
            new_children[child_name] = child
        else:
            if child_name is None:
                child_name = type(child).__name__

            found = False
            for i in range(1000):
                tmp_name = child_name + "_" + str(i)
                if tmp_name not in new_children:
                    found = True
                    new_children[tmp_name] = child
                    break

            if not found:
                raise OverflowError(
                    "Could not patch child name without name collision after 1000 modifications"
                )
        return self.construct_copy(children=new_children)

    def assign_children(self, new_children: Mapping[Hashable, ChildType]) -> Self:
        """Helper function to assign new children to this node without changing the child or data type of the tree

        Unlike calling `construct_copy()` directly, this will retain already existing children under this node if `new_children` does not overwrite all keys
        in this node

        Parameters
        ----------
        new_children : Mapping[Hashable, ChildType]
            The mapping of *additional* children to be appended to this node's list of children.

        Returns
        -------
        Self
            A copy of this node but with potentially more or different child nodes.
        """
        # TODO: FIXME: Implement
        all_children = dict(self._children)
        all_children.update(new_children)
        return self.construct_copy(children=all_children)

    def is_level(self, target_level: str | list[str]) -> bool:
        """Check whether we are at a certain level in the ShnitselDB structure

        Parameters
        ----------
        target_level : str | Iterable[str]
            Desired level(s) to check for and accept as the target level.

        Returns
        -------
        bool
            True if the current node is of the required level or one of the required levels
        """
        if isinstance(target_level, list):
            return self._level_name in target_level
        else:
            return self._level_name == target_level

    @overload
    def collect_data(
        self,
        with_path: Literal[True],
    ) -> Iterator[tuple[str, DataType]]: ...

    @overload
    def collect_data(self, with_path: Literal[False] = False) -> Iterator[DataType]: ...

    def collect_data(
        self, with_path: bool = False
    ) -> Iterator[DataType] | Iterator[tuple[str, DataType]]:
        """Function to retrieve all data entries in the tree underneath this node.

        Helpful for aggregating across all entries in a subtree without the need for
        full hierarchical information.

        Parameters
        ----------
        with_path : bool, default=False
            Flag to obtain an iterable over the pairs of paths and data instead.

        Yields
        ------
        Iterator[Iterable[DataType]]
            An iterator over all the data entries in this subtree.
        Iterator[tuple[str, DataType]]
            An iterator over all the data entries in this subtree paired with their paths in the tree.
        """

        if with_path:
            if self.has_data:
                yield (self.path, self.data)
            for child in self.children.values():
                if child is not None:
                    yield from child.collect_data(with_path=True)
        else:
            if self.has_data:
                yield self.data
            for child in self.children.values():
                if child is not None:
                    yield from child.collect_data()

    def apply_data_attributes(
        self, properties: dict
    ) -> "Self | TreeNode[Any, DataType] | None":
        """

        Parameters
        ----------
        properties : dict
            The attributes to set with their respective values.

        Returns
        -------
        Self | TreeNode[Any, DataType]
            The subtree after the update
        """

        props = {}

        for k, v in properties.items():
            if v is not None:
                props[k] = v

        def update_attrs(data: DataType, _props: dict) -> DataType:
            if hasattr(data, 'attrs'):
                getattr(data, 'attrs').update(props)
            return data

        return self.map_data(lambda x: update_attrs(x, props))

    def map_flat_group_data(
        self, map_func: Callable[[Iterable[DataType]], ResType | None]
    ) -> "Self|TreeNode[Any, ResType]":
        """Helper function to apply a mapping function to all flat group nodes.

        Will only apply the mapping function to nodes of type `DataGroup` and only those who have exclusively `DataLeaf` children.

        Parameters
        ----------
        map_func : Callable[[Iterable[DataType]], ResType  |  None]
            Function mapping the data in the flat groups to a new result type

        Returns
        -------
        Self | TreeNode[Any, ResType]
             A new subtree structure, which will hold leaves with ResType data underneath each mapped group.
        """
        from .data_group import DataGroup
        from .data_leaf import DataLeaf

        def extended_mapper(
            flat_group: TreeNode[Any, DataType],
        ) -> TreeNode[Any, ResType]:
            assert isinstance(flat_group, DataGroup)
            assert len(flat_group.subgroups) == 0
            child_data = {k: v.data for k, v in flat_group.subleaves.items()}
            # Actually perform the mapping over child data
            res = map_func(child_data.values())

            new_leaf = DataLeaf[ResType](name="reduced", data=res)
            new_group = flat_group.construct_copy(children={new_leaf.name: new_leaf})
            return new_group

        def filter_flat_groups(node: TreeNode[Any, DataType]) -> bool:
            return isinstance(node, DataGroup) and node.is_flat_group

        return self.map_filtered_nodes(
            filter_func=filter_flat_groups, map_func=extended_mapper
        )  # type: ignore

    def group_data_by_metadata(self) -> Self | None:
        """Helper function to allow for grouping of data within the tree by the metadata
        extracted from Trajectories.

        Should only be called on trees where `DataType=Trajectory` or `DataType=Frames` or subtypes thereof.
        Will fail due to an attribute error or yield an empty tree otherwise.

        Returns
        -------
        Self
            A tree where leaves are grouped to have similar metadata and only leaves with the same metadata are within the same gorup.
        """
        return self.group_children_by(
            key_func=_trajectory_key_func, group_leaves_only=True
        )

    @property
    def as_stacked(self) -> MultiSeriesStacked | DataType:
        return self.to_stacked()

    def to_stacked(
        self, only_direct_children: bool = False
    ) -> MultiSeriesStacked | DataType:
        """Stack the trajectories in a subtree into a multi-trajetctory dataset.

        The resulting dataset has a new `frame` dimension along which we can iterate through all individual frames of all trajectories.

        Parameters
        ----------
        only_direct_children : bool, optional
            Whether to only gather trajectories from direct children of this subtree.

        Returns
        -------
        MultiSeriesStacked
            The resulting multi-trajectory dataset stacked along a `frame` dimension
        DataType
            If it is an xarray.DataArray tree that we are concatenating.
        """
        from shnitsel.data.shnitsel_db_helpers import concat_subtree
        import xarray as xr

        res = concat_subtree(self, only_direct_children)
        if isinstance(res, xr.Dataset):
            return MultiSeriesStacked(res)
        return res

    @property
    def as_layered(self) -> MultiSeriesLayered:
        return self.to_layered()

    def to_layered(self, only_direct_children: bool = False) -> MultiSeriesLayered:
        """Lazer the trajectories in a subtree into a multi-trajectory dataset.

        The resulting dataset has a new `trajectorz` dimension along which we can iterate through all individual frames of all trajectories.

        Parameters
        ----------
        only_direct_children : bool, optional
            Whether to only gather trajectories from direct children of this subtree.

        Returns
        -------
        MultiSeriesLayered
            The resulting multi-trajectory dataset layered along a `trajectory` dimension
        """
        from shnitsel.data.shnitsel_db_helpers import layer_subtree

        # TODO: FIXME: Convert to appropriate return type
        return MultiSeriesLayered(layer_subtree(self, only_direct_children))

    def sel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance: int | float | Iterable[int | float] | None = None,
        drop: bool = False,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with each array indexed by tick labels
        along the specified dimension(s).

        In contrast to `Dataset.isel`, indexers for this method should use
        labels instead of integers.

        Under the hood, this method is powered by using pandas's powerful Index
        objects. This makes label based indexing essentially just as fast as
        using integer indexing.

        It also means this method uses pandas's (well documented) logic for
        indexing. This means you can use string shortcuts for datetime indexes
        (e.g., '2000-01' to select all values in January 2000). It also means
        that slices are treated as inclusive of both the start and stop values,
        unlike normal Python indexing.

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by scalars, slices or arrays of tick labels. For dimensions with
            multi-index, the indexer may also be a dict-like object with keys
            matching index level names.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method to use for inexact matches:

            * None (default): only exact matches
            * pad / ffill: propagate last valid index value forward
            * backfill / bfill: propagate next valid index value backward
            * nearest: use nearest valid index value
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
        drop : bool, optional
            If ``drop=True``, drop coordinates variables in `indexers` instead
            of making them scalar.
        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            A new Dataset with the same contents as this dataset, except each
            variable and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        :func:`Dataset.isel <Dataset.isel>`
        :func:`DataArray.sel <DataArray.sel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        """
        raise NotImplementedError(
            ".sel() not yet implemented for %s in hierarchical tree structures"
            % type(self)
        )
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "sel")
        query_results = map_index_queries(
            self, indexers=indexers, method=method, tolerance=tolerance
        )

        if drop:
            no_scalar_variables = {}
            for k, v in query_results.variables.items():
                if v.dims:
                    no_scalar_variables[k] = v
                elif k in self._coord_names:
                    query_results.drop_coords.append(k)
            query_results.variables = no_scalar_variables

        result = self.isel(indexers=query_results.dim_indexers, drop=drop)
        return result._overwrite_indexes(*query_results.as_tuple()[1:])

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        drop: bool = False,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new tree indexed along dimensions `compound`, `group` or `trajectory`
        and with data in leaves of the tree indexed along the remaining specified
        dimension(s) if the leaves support `.isel()` operations.

        Internally, it filters data with their own `.isel()` functions and performs
        some additional filtering specific to the tree structure

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by integers, slice objects or arrays.
            indexer can be a integer, slice, array-like or DataArray.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        drop : bool, default: False
            If ``drop=True``, drop coordinates variables indexed by integers
            instead of making them scalar.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Dataset:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : TreeNode[ChildType, DataType]
            A new tree with the same contents as this tree, except each
            data entry is indexed by the appropriate indexers and subtrees are filtered
            by the choices in tree-specific dimensions.
            The logic for selection on the leaf data entries is specific to the type of data in the leaf.

        Examples
        --------
        # TODO: FIXME: Provide better tree selection example.

        >>> tree = xr.Dataset(
        ...     {
        ...         "math_scores": (
        ...             ["student", "test"],
        ...             [[90, 85, 92], [78, 80, 85], [95, 92, 98]],
        ...         ),
        ...         "english_scores": (
        ...             ["student", "test"],
        ...             [[88, 90, 92], [75, 82, 79], [93, 96, 91]],
        ...         ),
        ...     },
        ...     coords={
        ...         "student": ["Alice", "Bob", "Charlie"],
        ...         "test": ["Test 1", "Test 2", "Test 3"],
        ...     },
        ... )

        # A specific element from the dataset is selected

        >>> dataset.isel(student=1, test=0)
        <xarray.Dataset> Size: 68B
        Dimensions:         ()
        Coordinates:
            student         <U7 28B 'Bob'
            test            <U6 24B 'Test 1'
        Data variables:
            math_scores     int64 8B 78
            english_scores  int64 8B 75

        # Indexing with a slice using isel

        >>> slice_of_data = dataset.isel(student=slice(0, 2), test=slice(0, 2))
        >>> slice_of_data
        <xarray.Dataset> Size: 168B
        Dimensions:         (student: 2, test: 2)
        Coordinates:
          * student         (student) <U7 56B 'Alice' 'Bob'
          * test            (test) <U6 48B 'Test 1' 'Test 2'
        Data variables:
            math_scores     (student, test) int64 32B 90 85 78 80
            english_scores  (student, test) int64 32B 88 90 75 82

        # Indexing using a sequence of keys.

        See Also
        --------
        :func:`Dataset.isel <Dataset.isel>`
        :func:`TreeNode.sel <TreeNode.sel>`

        """
        raise NotImplementedError(
            ".isel() not yet implemented for %s in hierarchical tree structures"
            % type(self)
        )
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")
        if any(is_fancy_indexer(idx) for idx in indexers.values()):
            return self._isel_fancy(indexers, drop=drop, missing_dims=missing_dims)

        # Much faster algorithm for when all indexers are ints, slices, one-dimensional
        # lists, or zero or one-dimensional np.ndarray's
        indexers = drop_dims_from_indexers(indexers, self.dims, missing_dims)

        variables = {}
        dims: dict[Hashable, int] = {}
        coord_names = self._coord_names.copy()

        indexes, index_variables = isel_indexes(self.xindexes, indexers)

        for name, var in self._variables.items():
            # preserve variable order
            if name in index_variables:
                var = index_variables[name]
            else:
                var_indexers = {k: v for k, v in indexers.items() if k in var.dims}
                if var_indexers:
                    var = var.isel(var_indexers)
                    if drop and var.ndim == 0 and name in coord_names:
                        coord_names.remove(name)
                        continue
            variables[name] = var
            dims.update(zip(var.dims, var.shape, strict=True))

        return self._construct_direct(
            variables=variables,
            coord_names=coord_names,
            dims=dims,
            attrs=self._attrs,
            indexes=indexes,
            encoding=self._encoding,
            close=self._close,
        )

    def __str__(self) -> str:
        """A basic representation of this node.

        Only contains rudimentary information about this node. Use `repr()` for a more extensive representation.

        Returns
        -------
        str
            A string representation with minimal information.
        """
        params = {}
        if self.has_data:
            params['data'] = str(self.data)
        if self._level_name is not None:
            params['level'] = str(self._level_name)
        if self._children:
            child_keys = list(str(x) for x in self._children.keys())
            params['children'] = f"{len(self._children)}: " + "; ".join(child_keys)

        return f"{type(self)} [{params}]"

    def __repr__(self) -> str:
        """A simple representation of the data and structure of this subtree.

        _extended_summary_

        Returns
        -------
        str
            A string representation with more extensive information than that returned by `__str__()`
        """
        params = {}
        if self.has_data:
            params['data'] = repr(self.data)
        if self._level_name is not None:
            params['level'] = str(self._level_name)
        if self._children:
            childrep = {k: repr(x) for k, x in self._children.items()}
            params['children'] = f"{len(self._children)}: " + repr(childrep)

        return f"{type(self)} [{params}]"

    def _repr_html_(self) -> str:
        """Obtain an html representation of this subtree.

        Currently generates a tabular representation of the subtree.

        Returns
        -------
        str
            A html string representing the data in this subtree.
        """
        # res = f"<h1>{type(self).__name__} (level: {self._level_name or 'unknown'})</h1>"
        # if self.has_data:
        #     res += "<br/>\r\n<h2>Data</h2>\r\n"
        #     sdata = self.data
        #     if hasattr(sdata, '_repr_html_') and False:
        #         res += f"<p>{sdata._repr_html_()}</p>\r\n"
        #     else:
        #         res += f"<p>{type(sdata).__qualname__}</p>\r\n"
        # if self._children:
        #     cres = "<h2>Children</h2>\r\n<table>\r\n"

        #     for k, c in self._children.items():
        #         if c is not None:
        #             cres += f"<tr><td>{repr(k)}</td><td>{c._repr_html_()}</td></tr>\r\n"
        #     cres += "</table>\r\n"
        #     res += cres
        # attrs = self.attrs
        # if attrs:
        #     ares = "<br/>\r\n<h2>Attributes</h2>\r\n<table>\r\n"

        #     for k, c in attrs.items():
        #         if c is not None:
        #             ares += f"<tr><td>{repr(k)}</td><td>{repr(c)}</td></tr>\r\n"
        #     ares += "</table>\r\n"
        #     res += ares
        # return "<div>" + res + "</div>"
        from .tree_vis import tree_repr

        # TODO: Consider options: https://github.com/etetoolkit/ete, https://treelib.readthedocs.io/en/latest/, https://plotly.com/python/tree-plots/
        return tree_repr(self)


def _trajectory_key_func(node: TreeNode) -> None | str | TrajectoryGroupingMetadata:
    """Helper function to extract trajectory metadata of leaf nodes for trees with
    appropriate data types.

    If applied to other nodes may yield a `None` key or just their `name` attribute as a `str`.

    Parameters
    ----------
    node : TreeNode
        The node to extract the `TrajectoryGroupingMetadata` metadata from.
        See `Trajectory.get_grouping_metadata()` for creation of the meta data
        instance.

    Returns
    -------
    None | str | TrajectoryGroupingMetadata
        The key to use for the grouping of this node.
    """
    from .data_leaf import DataLeaf
    from ..dataset_containers import Trajectory, Frames

    if isinstance(node, DataLeaf):
        if not node.has_data:
            # Do not group empty data
            return None
        else:
            # Get grouping metadata
            if isinstance(node.data, Trajectory) or isinstance(node.data, Frames):
                return node.data.get_grouping_metadata()
        # Don't attempt to group weird data types
        return None
    else:
        return node.name
