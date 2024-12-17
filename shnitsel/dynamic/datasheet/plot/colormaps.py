import matplotlib as mpl
import numpy as np

__all__ = ['magma_rw', 'custom_ylgnr']

_clmagma = mpl.colormaps["magma_r"](np.linspace(0, 1, 128))
_clmagma[:, 2] *= 1 / np.max(
    _clmagma[:, 2]
)  # more blue near zero, so white rather than yellow
magma_rw = mpl.colors.LinearSegmentedColormap.from_list('magma_rw', _clmagma)

custom_ylgnr = mpl.colors.LinearSegmentedColormap.from_list(
    'custom', mpl.colormaps['YlGn_r'](np.linspace(0, 0.75, 128))
)