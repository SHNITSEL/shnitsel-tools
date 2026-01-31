import builtins

# TODO FIXME: The following constant helps prevent Segfault on Python 3.13.2
# Is there a better way of checking whether the interpreter
# supports IPython import?
HAS_IPYTHON: bool = hasattr(builtins, '__IPYTHON__')

DATASET_ACCESSOR_NAME: str
DATAARRAY_ACCESSOR_NAME: str
DATASET_ACCESSOR_REGISTERED: bool = False
DATAARRAY_ACCESSOR_REGISTERED: bool = False