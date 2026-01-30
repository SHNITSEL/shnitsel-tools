"""This submodule contains functions used to interface with other packages and programs, especially RDKit."""

import logging
from typing import Literal

import numpy as np
from rdkit import Chem as rc
from rdkit.Chem import rdDepictor

from shnitsel._contracts import needs
from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.rd import set_atom_props, mol_to_numbered_smiles
from .core.typedefs import AtXYZ
from .units.conversion import convert_length
from .units.definitions import length
import xarray as xr


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
        logging.warning(
            "da.attrs['units'] is not set, the output will contain unconverted values"
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
        logging.warning(
            "da.attrs['units'] is not set, the output will contain unconverted values"
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
    atXYZ_frame: AtXYZ | xr.Dataset | ShnitselDataset,
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
    atXYZ_da_frame: xr.DataArray
    if isinstance(atXYZ_frame, xr.Dataset):
        atXYZ_da_frame = atXYZ_frame.atXYZ
    elif isinstance(atXYZ_frame, ShnitselDataset):
        atXYZ_da_frame = atXYZ_frame.positions
    elif isinstance(atXYZ_frame, xr.DataArray):
        atXYZ_da_frame = atXYZ_frame
    else:
        raise ValueError(
            "Unsupported type for `atXYZ_frame` parameter: %s" % type(atXYZ_frame)
        )

    # Make sure the unit is correct
    atXYZ_in_angstrom = convert_length(atXYZ_da_frame, to=length.Angstrom)
    mol = rc.rdmolfiles.MolFromXYZBlock(to_xyz(atXYZ_in_angstrom))
    rc.rdDetermineBonds.DetermineConnectivity(mol, useVdw=True, covFactor=covFactor)
    try:
        rc.rdDetermineBonds.DetermineBondOrders(mol, charge=(charge or 0))
    except ValueError:
        if charge is not None:
            raise
    if to2D:
        rdDepictor.Compute2DCoords(mol)  # type: ignore
    return set_atom_props(
        mol, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel
    )


def numbered_smiles_to_mol(smiles: str) -> rc.Mol:
    """Convert a numbered SMILES-string to a analogically-numbered Mol object

    Parameters
    ----------
    smiles : str
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
    mol = rc.RenumberAtoms(mol, map_new_to_old)
    return set_atom_props(mol, molAtomMapNumber=False)


def _most_stable_frame(atXYZ, obj: xr.Dataset | ShnitselDataset) -> xr.DataArray:
    """Find the frame, out of all the initial conditions,
    with the lowest ground-state energy;
    failing that, return the first frame in ``atXYZ``
    """
    leading_dim: str
    if isinstance(obj, DataSeries):
        leading_dim = obj.leading_dimension
    else:
        if 'frame' in atXYZ.sizes:
            leading_dim = 'frame'
        elif 'time' in atXYZ.sizes:
            leading_dim = 'time'
        else:
            logging.info("No leading dimension detected for atXYZ source")
            return atXYZ

    if (
        'energy' not in obj
        or 'state' not in obj
        or 'time' not in obj
        or (
            'trajid' not in obj.coords
            and 'trajid_' not in obj.coords
            and 'trajectory' not in obj.coords
            and 'atrajectory' not in obj.coords
        )
    ):
        return atXYZ.isel({leading_dim: 0})

    try:
        if 'trajectory' in obj.sizes or 'atrajectory' in obj.sizes:
            inicond_energy = obj.energy.isel(state=0).groupby('atrajectory').first()
            trajid = int(inicond_energy.idxmin().item())
            # TODO: FIXME: Why do we lack an index for `atrajectory`?
            return atXYZ.sel({'atrajectory': trajid}).isel(time=0)
        elif 'trajid' in obj.coords:
            inicond_energy = obj.energy.isel(state=0).groupby('trajid').first()
            trajid = int(inicond_energy.idxmin().item())
            return atXYZ.sel({'trajid': trajid, 'time': 0})
    except Exception as e:
        logging.debug(f"Failed detection of optimum frame for molecule: {e}")

    return atXYZ.isel({leading_dim: 0})


def construct_default_mol(
    obj: xr.Dataset | xr.DataArray | ShnitselDataset | rc.Mol,
    to2D: bool = True,
    charge: int | float | None = None,
    molAtomMapNumber: list[str] | Literal[True] | None = None,
    atomNote: list[str] | Literal[True] | None = None,
    atomLabel: list[str] | Literal[True] | None = None,
    silent_mode: bool = False,
) -> rc.Mol:
    """Try many ways to get a representative Mol object for an ensemble:

        1. Use the ``mol`` attr (of either obj or obj['atXYZ']) directly
        2. Feed the ``smiles_map`` attr (of either ``obj`` or ``obj['atXYZ']``) to
        :py:func:`shnitsel.bridges.default_mol`
        3. Take the geometry from the first frame of the molecule and the charge specified in the
        ``charge`` attr (charge=0 assumed if not specified) and feed these to
        :py:func:`shnitsel.bridges.to_mol`

    Parameters
    ----------
    obj : xr.Dataset | xr.DataArray | Trajectory | Frames | rc.Mol
        An 'atXYZ' xr.DataArray with molecular geometries
        or an xr.Dataset containing the above as one of its variables
        or an rc.Mol object that will just be returned.
    to2D : bool, optional
        Discard 3D information and generate 2D conformer (useful for displaying), by default True
    charge: int, float or None, optional
        Optional parameter to set the charge of the molecule if not present within the molecule data.
        If provided as an int, will be interpreted as number of elemental charges.
        Float will be converted to int and interpreted the same way.
        If not provided, will attempt to extract charge info from the xarray or Mol object and
        default to 0 charge if none can be found.
    molAtomMapNumber : list[str] | Literal[True], optional
        Set the ``molAtomMapNumber`` properties to values provided in a list,
        or (if ``True`` is passed) set the properties to the respective atom indices
    atomNote : list[str] | Literal[True], optional
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomNote`` properties
    atomLabel : list[str] | Literal[True], optional
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomLabel`` properties
    silent_mode: bool, optional
        Flag to disable logging outputs. Used for internal constructions of molecular structures.
        By default, status of mol construction is logged (i.e. ``silent_mode=False``).

    Returns
    -------
        An rdkit.Chem.Mol object

    Raises
    ------
    ValueError
        If the final approach fails

    Notes
    -----
    If this function uses an existing ``Mol`` object, it returns a copy.
    One consequence is that the decoration parameters
    ``molAtomMapNumber``, ``atomNote`` and ``atomLabel``
    do not affect the existing ``Mol`` object.
    """
    charge_int: int | None = None
    atXYZ: xr.DataArray

    def sap(mol):
        nonlocal molAtomMapNumber, atomNote, atomLabel
        return set_atom_props(
            mol,
            molAtomMapNumber=molAtomMapNumber,
            atomNote=atomNote,
            atomLabel=atomLabel,
        )

    if isinstance(obj, xr.Dataset):
        # TODO: FIXME: Make these internal attributes with double underscores so they don't get written out.
        if '__mol' in obj.attrs:
            return sap(rc.Mol(obj.attrs['__mol']))
        elif 'atXYZ' in obj:  # We have a frames Dataset
            atXYZ = _most_stable_frame(obj['atXYZ'], obj)
        else:
            raise ValueError("Not enough information to construct molecule from object")
    elif isinstance(obj, ShnitselDataset):
        if '__mol' in obj.attrs:
            return rc.Mol(obj.attrs['__mol'])
        atXYZ = _most_stable_frame(obj.positions, obj)
        if charge is None:
            charge = obj.charge
            if not silent_mode:
                logging.debug(f'{charge=}')
    else:
        if '__mol' in obj.attrs:
            return sap(rc.Mol(obj.attrs['__mol']))
        atXYZ = obj  # We have an atXYZ DataArray

    if charge is not None:
        if isinstance(charge, float):
            charge = int(np.round(charge))
        charge_int = charge
    if charge_int is None and 'charge' in obj.coords:
        charge_int = int(np.round(obj.charge.item()))
    if charge_int is None and 'charge' in obj.attrs:
        charge_int = int(obj.attrs.get('charge', 0))
    if charge_int is None and 'state_charges' in obj.coords:
        charge_int = int(obj.state_charges[0].item())
    if charge_int is None:
        if not silent_mode:
            logging.info("Assuming molecular charge as 0")
        charge_int = 0

    # TODO: FIXME: Make these internal attributes with double underscores so they don't get written out.
    if '__mol' in atXYZ.attrs:
        return sap(rc.Mol(obj.attrs['__mol']))
    elif 'smiles_map' in obj.attrs:
        if not silent_mode:
            logging.debug("default_mol: Using `obj.attrs['smiles_map']`")
        mol = numbered_smiles_to_mol(obj.attrs['smiles_map'])
    elif 'smiles_map' in atXYZ.attrs:
        return sap(numbered_smiles_to_mol(atXYZ.attrs['smiles_map']))

    if 'frame' in atXYZ.dims:
        if not silent_mode:
            logging.info("Picking first frame for molecule construction")
        atXYZ = atXYZ.isel(frame=0)
        if 'frame' in atXYZ.dims:
            atXYZ = atXYZ.squeeze('frame')
    if 'time' in atXYZ.dims:
        if not silent_mode:
            logging.info("Picking first time step for molecule construction")
        atXYZ = atXYZ.isel(time=0)
        if 'time' in atXYZ.dims:
            atXYZ = atXYZ.squeeze('time')

    try:
        if charge_int != 0 and not silent_mode:
            logging.info(f"Creating molecule with {charge_int=}")
        return sap(to_mol(atXYZ, charge=charge_int, to2D=to2D))
    except (KeyError, ValueError) as e:
        if not silent_mode:
            logging.error(e)
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


default_mol = construct_default_mol
