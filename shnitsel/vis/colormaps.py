import matplotlib as mpl
from matplotlib.colors import Colormap
import numpy as np

st_grey = '#2c3e50'
st_yellow = '#C4A000'  # (196/255, 160/255, 0/255)
st_violet = '#7E5273'
st_blue = '#2233AA'
st_orange = '#F67E00'
st_pink = "#E00D99"

__all__ = ['magma_rw', 'custom_ylgnr']

_clmagma: Colormap = mpl.colormaps["magma_r"](np.linspace(0, 1, 128))
_clmagma[:, 2] *= 1 / np.max(
    _clmagma[:, 2]
)  # more blue near zero, so white rather than yellow
magma_rw: Colormap = mpl.colors.LinearSegmentedColormap.from_list('magma_rw', _clmagma)

custom_ylgnr: Colormap = mpl.colors.LinearSegmentedColormap.from_list(
    'custom', mpl.colormaps['YlGn_r'](np.linspace(0, 0.75, 128))
)

_cmap_mult_degeneracy: Colormap = mpl.colormaps["Greys"]

default_singlet_state_colormap: Colormap = mpl.colormaps.get_cmap("winter")
default_doublet_state_colormap: Colormap = mpl.colormaps.get_cmap("cool")
default_triplet_state_colormap: Colormap = mpl.colormaps.get_cmap("autumn")

default_lower_singlet_colors = [
    '#000000',
    '#2c3e50',
    '#F6BE00',
    '#7515AD',
    # '#2233AA',
    #     st_yellow,
    #     st_grey,
    #     st_violet,
    # ]  #
    # '#4DAD15',
    # '#AD2915',
    # '#7515AD',
    # '#FF4D00',
]

default_lower_singlet_transition_colors = {
    (1, 2): '#F6BE00',
    (1, 3): '#2233AA',
    (2, 3): '#4DAD15',
}

default_lower_triplet_transition_colors = {
    (0, 1): '#AD2915',
    (0, 2): '#FF4D00',
    (1, 2): '#F67E00',
}


def get_default_singlet_state_colormap(num_singlets: int) -> list:
    """Get a list of per-state colors for this number of singlets.

    Parameters
    ----------
    num_singlets : int
            The number of maximum singlet states to consider

    Returns
    -------
    list
        The list of colors for each of these states
    """
    if num_singlets <= len(default_lower_singlet_colors):
        colors = [hex2rgb(x) for x in default_lower_singlet_colors][:num_singlets]
    else:
        rem_states = num_singlets - len(default_lower_singlet_colors)
        rem_colors = default_singlet_state_colormap(
            np.linspace(0, 0.85, num=rem_states)
        )
        colors = [hex2rgb(x) for x in default_lower_singlet_colors] + list(rem_colors)
    return colors


def get_default_doublet_state_colormap(
    num_doublets: int,
    degeneracy_groups: list[int] | None = None,
) -> list:
    """Get a list of per-state colors for this number of doublets.

    Parameters
    ----------
    num_doublets : int
            The number of doublets to generate colors for
    degeneracy_groups : dict[int, int], optional
            Degeneracy group indices to make sure different groups have the most distinct colors. If not set, all states will be considered different.s

    Returns
    -------
    list
        The list of colors for each of these states
    """
    # colors = list(default_doublet_state_colormap(np.linspace(0, 1.0, num=num_doublets)))
    return _get_default_degenerate_state_map(
        num_doublets, default_doublet_state_colormap, degeneracy_groups
    )


def get_default_triplet_state_colormap(
    num_triplets: int,
    degeneracy_groups: list[int] | None = None,
) -> list:
    """Get a list of per-state colors for this number of triplets.

    Parameters
    ----------
    num_triplets : int
            The number of triplets to generate colors for
    degeneracy_groups : dict[int, int], optional
            Degeneracy group indices to make sure different groups have the most distinct colors. If not set, all states will be considered different.

    Returns
    -------
    list
        The list of colors for each of these states
    """
    # colors = list(default_triplet_state_colormap(np.linspace(0, 1.0, num=num_triplets)))
    return _get_default_degenerate_state_map(
        num_triplets, default_triplet_state_colormap, degeneracy_groups
    )


# Coefficients and bias colors for small degrees of degeneracy
color_mix_info: dict[int, list[tuple[float, np.ndarray]]] = {
    1: [(0.0, np.array([1.0, 1.0, 1.0, 1.0]))],
    2: [(0.0, np.array([1.0, 1.0, 1.0, 1.0])), (0.3, np.array([1.0, 1.0, 1.0, 1.0]))],
    3: [
        (0.0, np.array([1.0, 1.0, 1.0, 1.0])),
        (0.3, np.array([0.0, 0.0, 0.0, 1.0])),
        (0.3, np.array([1.0, 1.0, 1.0, 1.0])),
    ],
}


def _get_default_degenerate_state_map(
    num_states: int, colormap: Colormap, degeneracy_groups: list[int] | None = None
) -> list:
    """Helper function to try and spread colors most strongly between different degeneracy groups.

    Parameters
    ----------
    num_states : int
            Number of states in total
    colormap : Colormap
            The colormap to use for the interstate colors
    degeneracy_groups : list[int] | None, optional
            The index of the degeneracy group for each state, only if there may be overlap. Defaults to None.

    Returns
    -------
    list
        The default mapping for degenerate states
    """

    if degeneracy_groups is not None:
        assert num_states == len(degeneracy_groups), (
            f"discrepancy between #states {num_states} and {len(degeneracy_groups)=}"
        )

        degeneracy_group_set = list(set(degeneracy_groups))
        degeneracy_group_set.sort()
        num_deg_groups = len(degeneracy_group_set)

        group_colorlist = list(colormap(np.linspace(0, 1.0, num=num_deg_groups)))

        degeneracy_group_states = {}
        colorlist_reordered = [None] * num_states

        for index, group in enumerate(degeneracy_groups):
            if group not in degeneracy_group_states:
                degeneracy_group_states[group] = []
            degeneracy_group_states[group].append(index)

        for group, g_color in zip(degeneracy_group_set, group_colorlist):
            group_states = degeneracy_group_states[group]
            group_states.sort()

            size_group = len(group_states)

            if size_group in color_mix_info:
                mix_data = color_mix_info[size_group]

                for index, (weight, offset) in zip(group_states, mix_data):
                    colorlist_reordered[index] = (g_color + weight * offset) / (
                        1.0 + weight
                    )
            else:
                for index, color_offset in zip(
                    group_states,
                    _cmap_mult_degeneracy(
                        colormap(np.linspace(0, 1.0, num=size_group))
                    ),
                ):
                    coefficient = 0.3
                    colorlist_reordered[index] = (
                        g_color + coefficient * color_offset
                    ) / (1 + coefficient)

        return colorlist_reordered
    else:
        colorlist = list(colormap(np.linspace(0, 1.0, num=num_states)))
        return colorlist


def get_default_state_colormap(
    num_states: int,
    multiplicity: int = 1,
    degeneracy_groups: list[int] | None = None,
) -> list[str]:
    """Get default state colormap for a number of states of a specific multiplicity

    Parameters
    ----------
    num_states : int
            Number of states in this multiplicity
    multiplicity : int, optional
            The multiplicity to get the colors for. Defaults to 1, i.e. Singlets.
    degeneracy_groups : list[int], optional
            Degeneracy group indices to make sure different groups have the most distinct colors. If not set, all states will be considered different.

    Returns
    -------
    list[str]
        The string representations of per-state colors
    """
    base = None
    if multiplicity == 1:
        base = get_default_singlet_state_colormap(num_states)
    elif multiplicity == 3:
        base = get_default_triplet_state_colormap(
            num_states, degeneracy_groups=degeneracy_groups
        )
    else:
        base = get_default_doublet_state_colormap(
            num_states, degeneracy_groups=degeneracy_groups
        )

    return [rgb2hex(x) for x in base]


def hex2rgb(hex_str: str) -> np.ndarray:
    return np.array(mpl.colors.to_rgb(hex_str))


def rgb2hex(rgb: np.ndarray) -> str:
    return mpl.colors.to_hex(rgb)


multiplicity_intra_bias = {
    1: np.array([1.0, 1.0, 1.0]),  # np.array([0.5, 0, 0]),
    2: np.array([0, 0.5, 0]),
    3: np.array([1.0, 1.0, 1.0]),  # np.array([0, 0, 0.5]),
}


def get_default_interstate_colormap_same_mult(
    multiplicity: int, colors: dict[int, str]
) -> dict[tuple[int, int], str]:
    """Function to generate a default inter-state colormap between states of the same multiplicity

    Parameters
    ----------
    multiplicity : int
            The state multiplicity in this set
    colors : dict[int, str]
            The colors do use as a basis for the states individually.

    Returns
    -------
    dict[tuple[int,int], str]
        Default interstate colors for states within the same multiplicity
    """

    mapped_colors = {k: hex2rgb(v) for k, v in colors.items()}
    min_index = np.min(list(colors.keys()))
    max_index = np.max(list(colors.keys()))

    index_span = max(1, np.abs(max_index - min_index))

    res_map = {}
    for k1 in mapped_colors:
        for k2 in mapped_colors:
            bias_coeff = (index_span - abs(k2 - k1)) / index_span
            total_intra = (
                multiplicity_intra_bias[multiplicity] * bias_coeff
                + mapped_colors[k1]
                + mapped_colors[k2]
            ) / (2.0 + bias_coeff)
            res_map[(k1, k2)] = rgb2hex(total_intra)
    if multiplicity == 1:
        for sc, color in default_lower_singlet_transition_colors.items():
            if sc in res_map:
                res_map[sc] = color
    if multiplicity == 3:
        for sc, color in default_lower_triplet_transition_colors.items():
            shifted_sc = (min_index + sc[0], min_index + sc[1])
            if shifted_sc in res_map:
                res_map[shifted_sc] = color

    return res_map


def get_default_interstate_colormap_inter_mult(
    colors_mult_1: dict[int, str],
    colors_mult_2: dict[int, str],
) -> dict[tuple[int, int], str]:
    """Function to generate a default inter-state colormap between states of the differents multiplicity

    Parameters
    ----------
    colors_mult_1 : dict[int, str]
            State color map of the first state multiplicity
    colors_mult_2 : dict[int, str]
            State color map of the second state multiplicity

    Returns
    -------
    dict[tuple[int,int], str]
        Resulting inter-state colormap
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
