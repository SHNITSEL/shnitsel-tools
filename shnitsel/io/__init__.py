

from shnitsel.core.postprocess import broaden_gauss as broaden_gauss
from shnitsel.core.parse import read_trajs as read_trajs
from shnitsel.core.ase import read_ase as read_ase

from .read import read
from .write import write

# Backward compatibility
from .parse_shnitsel import read_shnitsel_file as open_frames

__all__ = ['read', 'write', 'open_frames', 'read_trajs', 'read_ase']
