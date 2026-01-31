from functools import wraps
from typing import Callable, ParamSpec, TypeVar


Param = ParamSpec("Param")
RetType = TypeVar("RetType")


def internal():
    """Decorator to mark a function as for internal use only.

    Stability of the function interface is not guaranteed.
    May modify the docstring of the function to signify that it is indeed intended for internal use only.

    Returns
    -------
    The function with an internal annotation that its use should be internal.
    """

    def decorator(
        ds_func: Callable[Param, RetType],
    ) -> Callable[Param, RetType]:
        @wraps(ds_func)
        def wrap(*params, **kwargs):
            return ds_func(*params, **kwargs)

        return wrap

    return decorator


def API():
    """Decorator to mark a function as a part of the public API.

    Stability of the function interface will be maintained across minor versions unless stated explicitly
    May modify the docstring of the function to signify that it is indeed part of the public API

    Returns
    -------
    The function with an internal annotation that its use may be public.
    """

    def decorator(
        ds_func: Callable[Param, RetType],
    ) -> Callable[Param, RetType]:
        @wraps(ds_func)
        def wrap(*params, **kwargs):
            return ds_func(*params, **kwargs)

        return wrap

    return decorator


public = API
