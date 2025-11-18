from functools import wraps
from typing import Callable, Concatenate, Dict, Literal, ParamSpec, TypeVar
from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
import xarray as xr

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

Param = ParamSpec("Param")
RetType = TypeVar("RetType")
T = TypeVar("T", bound=xr.DataTree)


# NOTE: TODO: FIXME: This combines multiple things that should be split: Converting from Dataset to Tree input,
# Aggregation, unwrapping and registering with the ShnitselDB class.
# This should be seen as a work in progress and not yet finished.
def add_as_tree_method(
    cls: type[T] = ShnitselDB,
    aggregate_prior: Literal['all', 'compound', 'group'] | Callable | None = None,
    aggregate_post: Literal['all', 'compound', 'group'] | Callable | None = None,
    unwrap_single_result=False,
) -> Callable[
    [Callable[Concatenate[Trajectory, Param], RetType]],
    Callable[Concatenate[Trajectory, Param], RetType],
]:
    """Decorator to add a function to the Tree/Database version of a shnitsel Dataset.

    Automatically maps the function that applies to a dataset over the trajectories in the tree.
    Via additional arguments, you can specify, which kind of pre- and postprocessing should be performed on the database to support the function.

    Args:
        cls (type, optional):  The class to add the method to. Defaults to ShnitselDB.
        aggregate_prior (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Preprocessing method to apply. Option 1: specify the scope within which all trajectories should be aggregated, i.e.
            "all": use all the trajectories in the set as base of inputs to the function,
            "compound": use only the trajectories per compound group as base for input,
            "group": use only the trajectories within the same Group as base for input.
            altn). Defaults to None.
            Option 2: Provide an explicit pre-processing function to turn the tree structure into a different tree with potentially fewer datasets.
            Option 3: Perform no pre-processing by setting `None`.
            Defaults to `None.
        aggregate_post (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Same semantics as for `aggregate_prior`, but now the aggregation is applied to the tree after applying the wrapped function to all trajectories. Defaults to None.
        unwrap_single_result (bool, optional): Whether a single result should be returned as the unwrapped value (True) or contained in the tree structure. Defaults to False.

    Returns:
        Callable[
            [Callable[Concatenate[Trajectory, Param], RetType]],
            Callable[Concatenate[Trajectory, Param], RetType]
            ]: Returns a decorator that accepts a function with a trajectory as its first parameter and returns the function unchanged.
    """

    def decorator(
        ds_func: Callable[Concatenate[Trajectory, Param], RetType],
    ) -> Callable[Concatenate[Trajectory, Param], RetType]:
        # TODO: FIXME: Patch the annotations and documentation of the wrapper function compared to the original
        @wraps(ds_func)
        def wrapper(self, *args: Param.args, **kwargs: Param.kwargs):
            def simple_helper(ds: Trajectory) -> RetType:
                """We simply add this so that we can apply the function with the correct arguments to all trajectories.

                Args:
                    ds (Trajectory): The single trajectory to apply this method to

                Returns:
                    RetType: The result of the function `ds_func` applied to this dataset.s
                """
                return ds_func(ds, *args, **kwargs)

            # TODO: FIXME: Perform preprocessing.
            res = self.map_over_trajectories(simple_helper)
            # TODO: FIXME: Perform postprocessing.
            # TODO: FIXME: Perform unwrapping.
            return res

        setattr(cls, ds_func.__name__, wrapper)
        return ds_func

    return decorator


# This example code works but has issues with autocompletion:
#
# from shnitsel.data.shnitsel_db.db_function_decorator import add_as_tree_method
# from shnitsel.data.trajectory_format import Trajectory

# @add_as_tree_method()
# def db_ident(traj:Trajectory) ->Trajectory:
#     return traj

# obj = ShnitselDB()

# print("TESTEST")
# print(obj.__dict__)
# print(obj.db_ident)

# sys.exit(0)


# NOTE: This decorator is meant to allow the input of a datatree as first argument instead of a dataset
def dataset_to_tree_method(
    cls: type[T] = ShnitselDB,
    aggregate_prior: Literal['all', 'compound', 'group']
    | Callable[[T], T]
    | None = None,
    aggregate_post: Literal['all', 'compound', 'group']
    | Callable[[T], T]
    | None = None,
    unwrap_single_result: bool = False,
    parallel: bool = True,
) -> Callable[
    [Callable[Concatenate[Trajectory, Param], RetType]],
    Callable[Concatenate[Trajectory | T, Param], RetType],
]:
    """Decorator to add support for Tree/Database inputs when it originally only supports individual xr.Datasets.

    Automatically maps the function that applies to a dataset over the trajectories in the tree.
    Via additional arguments, you can specify, which kind of pre- and postprocessing should be performed on the database to support the function.

    Args:
        cls (type, optional): The class to add support for. Defaults to ShnitselDB.
        aggregate_prior (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Preprocessing method to apply. Option 1: specify the scope within which all trajectories should be aggregated, i.e.
            "all": use all the trajectories in the set as base of inputs to the function,
            "compound": use only the trajectories per compound group as base for input,
            "group": use only the trajectories within the same Group as base for input.
            altn). Defaults to None.
            Option 2: Provide an explicit pre-processing function to turn the tree structure into a different tree with potentially fewer datasets.
            Option 3: Perform no pre-processing by setting `None`.
            Defaults to `None.
        aggregate_post (Literal["all", "compound", "group"] | Callable[[cls], cls] | None, optional): Same semantics as for `aggregate_prior`, but now the aggregation is applied to the tree after applying the wrapped function to all trajectories. Defaults to None.
        unwrap_single_result (bool, optional): Whether a single result should be returned as the unwrapped value (True) or contained in the tree structure. Defaults to False.
        parallel (bool, optional): Whether application to different trajectories should occur in parallel. Defaults to True.

    Returns:
        Callable[
            [Callable[Concatenate[Trajectory, Param], RetType]],
            Callable[Concatenate[Trajectory|T, Param], RetType]
            ]: Returns a decorator that accepts a function with a trajectory as its first parameter and returns the function now supporting the cls type as a first parameter.
    """

    def decorator(
        ds_func: Callable[Concatenate[Trajectory, Param], RetType],
    ) -> Callable[Concatenate[Trajectory | T, Param], RetType]:
        # TODO: FIXME: Patch the annotations and documentation of the wrapper function compared to the original
        @wraps(ds_func)
        def wrapper(ds: Trajectory | T, *args: Param.args, **kwargs: Param.kwargs):
            def simple_helper(ds: Trajectory) -> RetType:
                """We simply add this so that we can apply the function with the correct arguments to all trajectories.

                Args:
                    ds (Trajectory): The single trajectory to apply this method to

                Returns:
                    RetType: The result of the function `ds_func` applied to this dataset.s
                """
                return ds_func(ds, *args, **kwargs)

            # TODO: FIXME: Perform preprocessing.
            res = ds.map_over_trajectories(simple_helper)
            # TODO: FIXME: Perform postprocessing.
            # TODO: FIXME: Perform unwrapping.
            return res

        return wrapper

    return decorator
