"""This submodule contains functions used to interface with other packages and programs, especially RDKit."""

from logging import warning
from typing import Literal

import numpy as np
from rdkit import Chem as rc

from shnitsel._contracts import needs
from shnitsel.rd import set_atom_props, mol_to_numbered_smiles
from .core.typedefs import AtXYZ
from .units.conversion import convert_length
from .units.definitions import length


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def to_xyz(da: AtXYZ, comment='#', units='angstrom') -> str:
    """Convert an xr.DataArray of molecular geometry to an XYZ string

    Parameters
    ----------
    da
        A molecular geometry -- should have dimensions 'atom' and 'direction'
    comment
        The comment line for the XYZ, by default '#'
    units
        The units to which to convert before creating the XYZ string

    Returns
    -------
        The XYZ data as a string

    Notes
    -----
        The units of the outputs will be the same as the array;
        consider converting to angstrom first, as most tools will expect this.
    """
    if 'units' not in da.attrs:
        warning(
            "da.attrs['units'] is not set, "
            "the output will contain unconverted values"
        )
    else:
        da = convert_length(da, to=units)
    atXYZ = da.transpose('atom', 'direction').values
    atNames = da.atNames.values
    sxyz = np.char.mod('% 23.15f', atXYZ)
    sxyz = np.squeeze(sxyz)
    sxyz = np.hstack((atNames.reshape(-1, 1), sxyz))
    sxyz = np.apply_along_axis(lambda row: ''.join(row), axis=1, arr=sxyz)
    return f'{len(sxyz):>12}\n  {comment}\n' + '\n'.join(sxyz)


@needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
def traj_to_xyz(traj_atXYZ: AtXYZ, units='angstrom') -> str:
    """Convert an entire trajectory's worth of geometries to an XYZ string

    Parameters
    ----------
    traj_atXYZ
        Molecular geometries -- should have dimensions 'atom' and 'direction'; should
        also be groupable by 'time' (i.e. either have a 'time' dimension or
        a 'time' coordinate)
    units
        The units to which to convert before creating the XYZ string

    Returns
    -------
        The XYZ data as a string, with time indicated in the comment line of each frame

    Notes
    -----
        The units of the outputs will be the same as the array;
        consider converting to angstrom first, as most tools will expect this.
    """
    if 'units' not in traj_atXYZ.attrs:
        warning(
            "da.attrs['units'] is not set, "
            "the output will contain unconverted values"
        )
    else:
        traj_atXYZ = convert_length(traj_atXYZ, to=units)

    atXYZ = traj_atXYZ.transpose(..., 'atom', 'direction').values
    if atXYZ.ndim == 2:
        atXYZ = atXYZ[None, :, :]
    assert len(atXYZ.shape) == 3
    atNames = traj_atXYZ.atNames.values
    sxyz = np.strings.mod('% 13.9f', atXYZ)
    sxyz = atNames[None, :] + sxyz[:, :, 0] + sxyz[:, :, 1] + sxyz[:, :, 2]
    atom_lines = np.broadcast_to([str(traj_atXYZ.sizes['atom'])], (sxyz.shape[0], 1))
    if 'time' in traj_atXYZ.coords:
        time_values = np.atleast_1d(traj_atXYZ.coords['time'])
        comment_lines = np.strings.mod('# t=%.2f', time_values)[:, None]
    else:
        comment_lines = np.broadcast_to([''], (sxyz.shape[0], 1))
    return '\n'.join(np.concat([atom_lines, comment_lines, sxyz], 1).ravel())


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def to_mol(
    atXYZ_frame: AtXYZ,
    charge: int | None = None,
    covFactor: float = 1.2,
    to2D: bool = True,
    molAtomMapNumber: list | Literal[True] | None = None,
    atomNote: list | Literal[True] | None = None,
    atomLabel: list | Literal[True] | None = None,
) -> rc.Mol:
    """Convert a single frame's geometry to an RDKit Mol object

    Parameters
    ----------
    atXYZ_frame
        The ``xr.DataArray`` object to be converted; must have 'atom' and 'direction' dims,
        must not have 'frame' dim.
    charge
        Charge of the molecule, used by RDKit to determine bond orders; if ``None`` (the default),
        this function will try ``charge=0`` and leave the bond orders undetermined if that causes
        an error; otherwise failure to determine bond order will raise an error.
    covFactor
        Scales the distance at which atoms are considered bonded, by default 1.2
    to2D
        Discard 3D information and generate 2D conformer (useful for displaying), by default True
    molAtomMapNumber
        Set the ``molAtomMapNumber`` properties to values provided in a list,
        or (if ``True`` is passed) set the properties to the respective atom indices
    atomNote
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomNote`` properties
    atomLabel
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomLabel`` properties

    Returns
    -------
        An RDKit Mol object

    Raises
    ------
    ValueError
        If ``charge`` is not ``None`` and bond order determination fails
    """
    # Make sure the unit is correct
    atXYZ_in_angstrom = convert_length(atXYZ_frame, to=length.Angstrom)
    mol = rc.rdmolfiles.MolFromXYZBlock(to_xyz(atXYZ_in_angstrom))
    rc.rdDetermineBonds.DetermineConnectivity(mol, useVdw=True, covFactor=covFactor)
    try:
        rc.rdDetermineBonds.DetermineBondOrders(mol, charge=(charge or 0))
    except ValueError as err:
        if charge is not None:
            raise err
    if to2D:
        rc.rdDepictor.Compute2DCoords(mol)  # type: ignore
    return set_atom_props(
        mol, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel
    )


def numbered_smiles_to_mol(smiles: str) -> rc.Mol:
    """Convert a numbered SMILES-string to a analogically-numbered Mol object

    Parameters
    ----------
    smiles
        A SMILES string in which each atom is associated with a mapping index,
        e.g. '[H:3][C:1]#[C:0][H:2]'

    Returns
    -------
        An :py:func:`rdkit.Chem.Mol` object with atom indices numbered according
        to the indices from the SMILES-string
    """
    mol = rc.MolFromSmiles(smiles, sanitize=False)  # sanitizing would strip hydrogens
    map_new_to_old = [-1 for i in range(mol.GetNumAtoms())]
    for atom in mol.GetAtoms():
        # Renumbering with e.g. [3, 2, 0, 1] means atom 3 gets new index 0, not vice-versa!
        map_new_to_old[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    return rc.RenumberAtoms(mol, map_new_to_old)


def default_mol(obj) -> rc.Mol:
    """Try many ways to get a representative Mol object for an ensemble:

        1. Use the ``mol`` attr (of either obj or obj['atXYZ']) directly
        2. Feed the ``smiles_map`` attr (of either ``obj`` or ``obj['atXYZ']``) to
        :py:func:`shnitsel.bridges.default_mol`
        3. Take the geometry from the first frame of the molecule and the charge specified in the
        ``charge`` attr (charge=0 assumed if not specified) and feed these to
        :py:func:`shnitsel.bridges.to_mol`

    Parameters
    ----------
    obj
        An 'atXYZ' xr.DataArray with molecular geometries
        or an xr.Dataset containing the above as one of its variables

    Returns
    -------
        An rdkit.Chem.Mol object

    Raises
    ------
    ValueError
        If the final approach fails
    """
    if 'atXYZ' in obj:  # We have a frames Dataset
        atXYZ = obj['atXYZ']
    else:
        atXYZ = obj  # We have an atXYZ DataArray

    if 'mol' in obj.attrs:
        return rc.Mol(obj.attrs['mol'])
    elif 'mol' in atXYZ.attrs:
        return rc.Mol(obj.attrs['mol'])
    elif 'smiles_map' in obj.attrs:
        return numbered_smiles_to_mol(obj.attrs['smiles_map'])
    elif 'smiles_map' in atXYZ.attrs:
        return numbered_smiles_to_mol(atXYZ.attrs['smiles_map'])

    try:
        charge = obj.attrs.get('charge', 0)
        return to_mol(atXYZ.isel(frame=0), charge=charge)
    except (KeyError, ValueError):
        raise ValueError(
            "Failed to get default mol, please set a smiles map. "
            "For example, if the compound has charge c and frame i "
            "contains a representative geometry, use "
            "frames.attrs['smiles_map'] = frames.atXYZ.isel(frame=i).st.get_smiles_map(charge=c)"
        )


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def smiles_map(atXYZ_frame, charge=0, covFactor=1.5) -> str:
    """Convert a geometry to a SMILES-string, retaining atom order

    Parameters
    ----------
    atXYZ_frame
        An xr.DataArray of molecular geometry
    charge, optional
        The charge of the molcule, by default 0
    covFactor, optional
        Scales the distance at which atoms are considered bonded, by default 1.5

    Returns
    -------
        A SMILES-string in which the mapping number indicates the order in which the
        atoms appeared in the input matrix, e.g. '[H:3][C:1]#[C:0][H:2]'
    """
    mol = to_mol(atXYZ_frame, charge=charge, covFactor=covFactor, to2D=True)
    return mol_to_numbered_smiles(mol)
