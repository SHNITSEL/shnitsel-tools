# from functools import wraps
# import logging
# from typing import Any, Callable, Concatenate, Literal, ParamSpec, TypeVar
# from shnitsel.data.shnitsel_db_helpers import concat_subtree, layer_subtree
# from shnitsel.data.tree import ShnitselDB, TreeNode
# from shnitsel.data.shnitsel_db_helpers import (
#     aggregate_xr_over_levels,
#     unwrap_single_entry_in_tree,
# )
# from shnitsel.data.dataset_containers import Trajectory, Frames
# import xarray as xr

# How functools updates a wrapper to look like the original.
# WRAPPER_ASSIGNMENTS = ('__module__', '__name__', '__qualname__', '__doc__',
#                        '__annotations__', '__type_params__')
# WRAPPER_UPDATES = ('__dict__',)
# def update_wrapper(wrapper,
#                    wrapped,
#                    assigned = WRAPPER_ASSIGNMENTS,
#                    updated = WRAPPER_UPDATES):
#     """Update a wrapper function to look like the wrapped function

#        wrapper is the function to be updated
#        wrapped is the original function
#        assigned is a tuple naming the attributes assigned directly
#        from the wrapped function to the wrapper function (defaults to
#        functools.WRAPPER_ASSIGNMENTS)
#        updated is a tuple naming the attributes of the wrapper that
#        are updated with the corresponding attribute from the wrapped
#        function (defaults to functools.WRAPPER_UPDATES)
#     """
#     for attr in assigned:
#         try:
#             value = getattr(wrapped, attr)
#         except AttributeError:
#             pass
#         else:
#             setattr(wrapper, attr, value)
#     for attr in updated:
#         getattr(wrapper, attr).update(getattr(wrapped, attr, {}))
#     # Issue #17482: set __wrapped__ last so we don't inadvertently copy it
#     # from the wrapped function when updating __dict__
#     wrapper.__wrapped__ = wrapped
#     # Return the wrapper so this can be used as a decorator via partial()
#     return wrapper

# Code the below decorator is loosely based on
# def add_method(cls):
#     def decorator(func):
#         @wraps(func)
#         def wrapper(*args, **kwargs):
#             return func(*args, **kwargs)
#         setattr(cls, func.__name__, wrapper)
#         return func
#     return decorator

# Param = ParamSpec("Param")
# RetType = TypeVar("RetType")
# DataType = TypeVar("DataType")
# T = TypeVar("T", bound=TreeNode)


# def add_as_tree_method(
#     cls: type[T] = ShnitselDB[Trajectory],
# ) -> Callable[
#     [Callable[Concatenate[T, Param], RetType]],
#     Callable[Concatenate[T, Param], RetType],
# ]:
#     """Decorator to add a function to the Tree/Database version of a shnitsel Dataset.

#     If the function does not accept the tree-type argument as first argument, it can be transformed to support it with the decorator `dataset_to_tree_method`
#     Parameters
#     ----------
#         cls (type, optional):  The class to add the method to. Defaults to ShnitselDB.

#     Returns
#     -------
#         Callable[
#             [Callable[Concatenate[T, Param], RetType]],
#             Callable[Concatenate[T, Param], RetType]
#             ]: Returns the function unchanged but adds it to the `cls` class as an instance method.
#     """

#     def decorator(
#         ds_func: Callable[Concatenate[T, Param], RetType],
#     ) -> Callable[Concatenate[T, Param], RetType]:
#         # TODO: FIXME: Make sure that this shows up in autocompletion and documentation
#         setattr(cls, ds_func.__name__, ds_func)
#         return ds_func

#     return decorator


# # This example code works but has issues with autocompletion:
# #
# # from shnitsel.data.shnitsel_db.db_function_decorator import add_as_tree_method
# # from shnitsel.data.trajectory_format import Trajectory

# # @add_as_tree_method()
# # def db_ident(traj:Trajectory) ->Trajectory:
# #     return traj

# # obj = ShnitselDB()

# # print("TESTEST")
# # print(obj.__dict__)
# # print(obj.db_ident)

# # sys.exit(0)



# # TODO: FIXME: This currently does not work. Typing of the parameters for trees and callbacks are off.
# def dataset_to_tree_method(
#     cls: type[T] = ShnitselDB[Trajectory],
#     aggregate_pre: Literal['all', 'compound', 'group'] | Callable[[T], T] | None = None,
#     aggregate_method_pre: Literal['concat', 'layer', 'list'] | None = None,
#     aggregate_post: Callable[[T], T] | None = None,
#     map_result_as_dict: bool = False,
#     unwrap_single_result: bool = False,
#     parallel: bool = True,
# ) -> Callable[
#     [Callable[Concatenate[Trajectory, Param], RetType]],
#     Callable[Concatenate[Trajectory | T, Param], RetType | T | dict],
# ]:
#     """Decorator to add support for Tree/Database inputs when it originally only supports individual xr.Datasets.

#     Automatically maps the function that applies to a dataset over the trajectories in the tree.
#     Via additional arguments, you can specify, which kind of pre- and postprocessing should be performed on the database to support the function.

#     Parameters
#     ----------
#         cls (type, optional): The class to add support for. Defaults to ShnitselDB.
#         aggregate_prior (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Preprocessing method to apply. Option 1: specify the scope within which all trajectories should be aggregated, i.e.
#             "all": use all the trajectories in the set as base of inputs to the function,
#             "compound": use only the trajectories per compound group as base for input,
#             "group": use only the trajectories within the same Group as base for input.
#             altn). Defaults to None.
#             Option 2: Provide an explicit pre-processing function to turn the tree structure into a different tree with potentially fewer datasets.
#             Option 3: Perform no pre-processing by setting `None`.
#             Defaults to `None.
#         aggregate_post (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Same semantics as for `aggregate_prior`, but now the aggregation is applied to the tree after applying the wrapped function to all trajectories. Defaults to None.
#         unwrap_single_result (bool, optional): Whether a single result should be returned as the unwrapped value (True) or contained in the tree structure. Defaults to False.
#         parallel (bool, optional): Whether application to different trajectories should occur in parallel. Defaults to True.

#     Returns
#     -------
#         Callable[
#             [Callable[Concatenate[Trajectory, Param], RetType]],
#             Callable[Concatenate[Trajectory|T, Param], RetType]
#             ]: Returns a decorator that accepts a function with a trajectory as its first parameter and returns the function now supporting the cls type as a first parameter.
#     """

#     def decorator(
#         ds_func: Callable[Concatenate[xr.Dataset, Param], RetType],
#     ) -> Callable[Concatenate[xr.Dataset | T, Param], RetType | T]:
#         # TODO: FIXME: Patch the annotations and documentation of the wrapper function compared to the original
#         @wraps(ds_func)
#         def wrapper(ds: Trajectory | T, *args: Param.args, **kwargs: Param.kwargs):
#             if isinstance(ds, cls):
#                 tree = ds

#                 def simple_helper(ds: xr.Dataset) -> RetType:
#                     """We simply add this so that we can apply the function with the correct arguments to all trajectories.

#                     Args:
#                         ds (Trajectory): The single trajectory to apply this method to

#                     Returns:
#                         RetType: The result of the function `ds_func` applied to this dataset.s
#                     """
#                     return ds_func(ds, *args, **kwargs)

#                 if aggregate_pre is not None:
#                     # Apply pre-processing
#                     if callable(aggregate_pre):
#                         tree = aggregate_pre(tree)
#                     else:
#                         if aggregate_method_pre == "concat":
#                             aggr_func = concat_subtree
#                         elif aggregate_method_pre == "layer":
#                             aggr_func = layer_subtree
#                         elif aggregate_method_pre == "list":
#                             raise NotImplementedError(
#                                 f"Aggregation method {aggregate_method_pre} requested at DB level {aggregate_pre} has not been implemented."
#                             )
#                             # aggr_func = list_subtree
#                         else:
#                             logging.error(
#                                 f"Unknown aggregation method {aggregate_method_pre} requested at DB level {aggregate_pre}"
#                             )
#                             raise ValueError(
#                                 f"Unknown aggregation method {aggregate_method_pre} requested at DB level {aggregate_pre}"
#                             )
#                         tree_res = aggregate_xr_over_levels(
#                             tree, aggr_func, aggregate_pre
#                         )
#                         if tree_res is None:
#                             raise ValueError(
#                                 f"Aggregation method {aggregate_method_pre} at DB level {aggregate_pre} yielded no result. Make sure that the tree is not empty."
#                             )
#                         else:
#                             tree = tree_res
#                 tree: T
#                 # TODO: FIXME: The return type of the function (Trajectory) does not match RetType. We may want to restrict RetType
#                 # Perform preprocessing.
#                 res = tree.map_data(
#                     simple_helper, parallel=parallel, result_as_dict=map_result_as_dict
#                 )
#                 # T
#                 # assert isinstance(res, cls)

#                 # Apply postprocessing
#                 if aggregate_post is not None:
#                     res = aggregate_post(res)

#                 # Unwrap result if requested
#                 if not map_result_as_dict and unwrap_single_result:
#                     res: Trajectory | T = unwrap_single_entry_in_tree(res)

#                 return res
#             else:
#                 assert isinstance(ds, Trajectory)
#                 return ds_func(ds, *args, **kwargs)

#         return wrapper

#     return decorator
