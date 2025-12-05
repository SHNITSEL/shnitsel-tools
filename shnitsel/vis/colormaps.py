import matplotlib as mpl
from matplotlib.colors import Colormap
import numpy as np

st_grey = '#2c3e50'
st_yellow = '#C4A000'  # (196/255, 160/255, 0/255)
st_violet = '#7E5273'

__all__ = ['magma_rw', 'custom_ylgnr']

_clmagma: Colormap = mpl.colormaps["magma_r"](np.linspace(0, 1, 128))
_clmagma[:, 2] *= 1 / np.max(
    _clmagma[:, 2]
)  # more blue near zero, so white rather than yellow
magma_rw: Colormap = mpl.colors.LinearSegmentedColormap.from_list('magma_rw', _clmagma)

custom_ylgnr: Colormap = mpl.colors.LinearSegmentedColormap.from_list(
    'custom', mpl.colormaps['YlGn_r'](np.linspace(0, 0.75, 128))
)

default_singlet_state_colormap: Colormap = mpl.colormaps.get_cmap("winter")
default_doublet_state_colormap: Colormap = mpl.colormaps.get_cmap("cool")
default_triplet_state_colormap: Colormap = mpl.colormaps.get_cmap("autumn")

default_lower_singlet_colors = [
    '#000000'
]  # ,'#4DAD15', '#AD2915', '#7515AD', '#FF4D00']


def get_default_singlet_state_colormap(num_singlets: int) -> list:
    """Get a list of per-state colors for this number of singlets.

    Args:
        num_singlets (int): The number of maximum singlet states to consider

    Returns:
        list: The list of colors for each of these states
    """
    if num_singlets <= len(default_lower_singlet_colors):
        colors = [hex2rgb(x) for x in default_lower_singlet_colors][:num_singlets]
    else:
        rem_states = num_singlets - len(default_lower_singlet_colors)
        rem_colors = default_singlet_state_colormap(np.linspace(0, 1.0, num=rem_states))
        colors = [hex2rgb(x) for x in default_lower_singlet_colors] + list(rem_colors)
    return colors


def get_default_doublet_state_colormap(num_doublets: int) -> list:
    """Get a list of per-state colors for this number of doublets.

    Args:
        num_doublets (int): The number of doublets to generate colors for

    Returns:
        list: The list of colors for each of these states
    """
    colors = list(default_doublet_state_colormap(np.linspace(0, 1.0, num=num_doublets)))
    return colors


def get_default_triplet_state_colormap(num_triplets: int) -> list:
    """Get a list of per-state colors for this number of triplets.

    Args:
        num_triplets (int): The number of triplets to generate colors for

    Returns:
        list: The list of colors for each of these states
    """
    colors = list(default_triplet_state_colormap(np.linspace(0, 1.0, num=num_triplets)))
    return colors


def get_default_state_colormap(num_states: int, multiplicity: int = 1) -> list[str]:
    """Get default state colormap for a number of states of a specific multiplicity

    Args:
        num_states (int): Number of states in this multiplicity
        multiplicity (int, optional): The multiplicity to get the colors for. Defaults to 1, i.e. Singlets.

    Returns:
        list[str]: The string representations of per-state colors
    """
    base = None
    if multiplicity == 1:
        base = get_default_singlet_state_colormap(num_states)
    elif multiplicity == 3:
        base = get_default_triplet_state_colormap(num_states)
    else:
        base = get_default_doublet_state_colormap(num_states)

    return [rgb2hex(x) for x in base]


def hex2rgb(hex_str: str) -> np.ndarray:
    return np.array(mpl.colors.to_rgb(hex_str))


def rgb2hex(rgb: np.ndarray) -> str:
    return mpl.colors.to_hex(rgb)


multiplicity_intra_bias = {
    1: np.array([0.5, 0, 0]),
    2: np.array([0, 0.5, 0]),
    3: np.array([0, 0, 0.5]),
}


def get_default_interstate_colormap_same_mult(
    multiplicity: int, colors: dict[int, str]
) -> dict[tuple[int, int], str]:
    """Function to generate a default inter-state colormap between states of the same multiplicity

    Args:
        multiplicity (int): _description_
        colors (dict[int, str]): _description_

    Returns:
        dict[tuple[int,int], str]: _description_
    """

    mapped_colors = {k: hex2rgb(v) for k, v in colors.items()}

    res_map = {}
    for k1 in mapped_colors:
        for k2 in mapped_colors:
            total_intra = (
                multiplicity_intra_bias[multiplicity]
                + mapped_colors[k1]
                + mapped_colors[k2]
            ) / 3.0
            res_map[(k1, k2)] = rgb2hex(total_intra)
    return res_map


def get_default_interstate_colormap_inter_mult(
    colors_mult_1: dict[int, str],
    colors_mult_2: dict[int, str],
) -> dict[tuple[int, int], str]:
    """Function to generate a default inter-state colormap between states of the differents multiplicity

    Args:
        colors_mult_1 (dict[int, str]): State color map of the first state multiplicity
        colors_mult_2 (dict[int, str]): State color map of the second state multiplicity

    Returns:
        dict[tuple[int,int], str]: Resulting inter-state colormap
    """

    mapped_colors_1 = {k: hex2rgb(v) for k, v in colors_mult_1.items()}
    mapped_colors_2 = {k: hex2rgb(v) for k, v in colors_mult_2.items()}

    res_map = {}
    for k1 in mapped_colors_1:
        for k2 in mapped_colors_2:
            inter_color = rgb2hex((mapped_colors_1[k1] + mapped_colors_2[k2]) / 2.0)
            res_map[(k1, k2)] = inter_color
            res_map[(k2, k1)] = inter_color
    return res_map
