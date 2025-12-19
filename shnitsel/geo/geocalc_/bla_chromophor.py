
@needs(dims={'atom', 'direction'})
def get_bla_chromophor(
    atXYZ: xr.DataArray,
    matches_or_mol: dict | Mol | None = None,
    mol: Mol | None = None,
    ang: Literal[False, 'deg', 'rad'] = False,
) -> xr.DataArray:
    """Identify triples of bonded atoms (using RDKit) and calculate every bond angle for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    matches_or_mol, optional
        A list containing information for each internal coordinate to be calculated.
        It may be convenient to use :py:func:`shnitsel.geo.geomatch.flag_angles`
        to create a dictionary in the correct format, and then customize it.
        Alternatively, you may supply an RDKit ``Mol`` object, which is passed to
        :py:func:`shnitsel.geo.geomatch.flag_bats`.
        If this argument is omitted altogether, an RDKit ``Mol`` object is generated
        using :py:func:`shnitsel.bridges.default_mol` and used as above.

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond angles with dimensions `frame` and `angle`.

    Raises
    ------
    UserWarning
        If both `matches` and `mol` are specified.
    """
    matches = _check_matches(matches_or_mol, atXYZ, fn=flag_bla_chromophor)['bonds']

    _, atom_idxs, bond_idxs, bond_types, fragment_objs = zip(*matches)

    assert all(len(x) == 2 for x in atom_idxs)
    assert all(len(x) == 1 for x in bond_idxs)
    assert all(len(x) == 1 for x in bond_types)

    single_idxs = []
    double_idxs = []
    for _, idxs, _, btype, _ in matches:
        if btype[0] == 1:
            single_idxs.append(idxs)
        else:
            double_idxs.append(idxs)

    sr0, sr1 = _positions(atXYZ, single_idxs)
    dr0, dr1 = _positions(atXYZ, double_idxs)

    data = dnorm(sr0 - sr1).mean('descriptor') - dnorm(dr0 - dr1).mean('descriptor')
    data = data.expand_dims('descriptor')

    # joined_atom_idxs = tuple(reduce(lambda x, y: set(x).union(y), atom_idxs))
    joined_atom_idxs = tuple(set(sum(bond_idxs, ())))
    format_str = 'BLA$_{' + ','.join(['%d'] * len(joined_atom_idxs)) + '}$'
    return _assign_descriptor_coords(
        data,
        # *std_args,
        [joined_atom_idxs],
        [sum(bond_idxs, ())],
        [sum(bond_types, ())],
        [rc.Mol()],
        format_str,
        per_atom_coords=False,
    )
