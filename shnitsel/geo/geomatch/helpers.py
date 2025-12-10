
def __level_to_geo(flag_level: int) -> str:
    """Helper function to map a flagging level in integer values to a flagging level as a geo feature string.

    Args:
        flag_level (int): The flagging level

    Returns:
        str: The string describing the geometry feature level, i.e. `bonds`, `angles` or `dihedrals`.
    """
    level_to_prop = {
        0: 'bonds',
        1: 'bonds',
        2: 'angles',
        3: 'dihedrals',
        4: 'dihedrals',
    }

    if flag_level < 0:
        flag_level = 0
    elif flag_level > 4:
        flag_level = 4

    return level_to_prop[flag_level]
