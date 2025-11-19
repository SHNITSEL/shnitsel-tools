import os
from typing import Any, Collection, Literal

from ase import Atoms
from ase.db import connect
import numpy as np
import xarray as xr

from shnitsel._contracts import needs
from shnitsel.data.trajectory_format import Trajectory


def _prepare_for_write_schnetpack(
    traj: Trajectory, leading_dim_name: Literal['frame', 'time']
) -> Trajectory:
    """Helper function to perform some preprocessing on the dataset before writing to a SchnetPack compatible database.

    Combines the dipole variables into one entry.

    Args:
        traj (Trajectory): The Dataset to transform into a SchnetPack conforming format.
        leading_dim_name (Literal['frame', 'time']): The name of the leading dimension identifying different frames within the dataset. Depending on the setup, this should be 'frame' or 'time'.

    Returns:
        Trajectory: The transformed dataset
    """
    # Recombine permanent and transition dipoles, as schnetpack expects
    dipoles: np.ndarray | xr.DataArray | None = None
    traj = traj.copy(deep=False)
    dip_attributes = {}
    if 'dipoles' in traj:
        dipoles = traj['dipoles']
        dip_attributes = traj['dipoles'].attrs
    elif 'dip_perm' in traj and 'dip_trans' in traj:
        dip_perm = (
            traj['dip_perm'].transpose(leading_dim_name, 'state', 'direction').data
        )
        dip_trans = (
            traj['dip_trans'].transpose(leading_dim_name, 'statecomb', 'direction').data
        )
        dipoles = np.concat((dip_perm, dip_trans.data), axis=1)
        dip_attributes = dip_perm.attrs

        del traj['dip_perm'], traj['dip_trans']
    elif 'dip_perm' in traj:
        dipoles = traj['dip_perm']
        dip_attributes = dipoles.attrs
        del traj['dip_perm']
    elif 'dip_trans' in traj:
        dipoles = traj['dip_trans']
        dip_attributes = dipoles.attrs
        del traj['dip_trans']

    if dipoles is not None:
        # Change some attributes before assigning
        dip_attributes["long_name"] = "Combined dipole entry (dip_perm and dip_trans)"
        dip_attributes["description"] = (
            "Combined dipole moment containing both permanent and transitional dipole information (if available)"
        )

        traj['dipoles'] = (
            [leading_dim_name, 'state_or_statecomb', 'direction'],
            dipoles,
            dip_attributes,
        )

    return traj


def _collect_metadata(traj: Trajectory) -> dict[str, Any]:
    """Helper function to generate the SPaiNN Metadata dict from a Trajectory struct.

    Extracts info from attributes and variables to set up the dict.

    Args:
        traj (Trajectory): The Dataset to extract the metadata from.

    Returns:
        dict[str, Any]: The resulting metadata dictionary.
    """
    # Define metadata information (dictionary)
    metadata: dict[str, Any] = {}

    if "trajectory_input_path" in traj.attrs:
        metadata['info'] = traj.attrs["trajectory_input_path"]

    if "est_level" in traj.attrs:
        metadata['ReferenceMethod'] = traj.attrs["est_level"]
        # (
        #     'SA3-CASSCF(2,2)'  # state-average CASSCF with 2 electrons in 2 orbitals
        # )

    metadata['_distance_unit'] = traj["atXYZ"].attrs["unit"]

    metadata['_property_unit_dict'] = {
        'energy': traj["energy"].attrs["unit"]
        if "energy" in traj
        else "1",  # 'Hartree',
        'forces': traj["forces"].attrs["unit"]
        if "forces" in traj
        else "1",  # 'Hartree/Bohr',
        'nacs': traj["nacs"].attrs["unit"]
        if "nacs" in traj
        else "1",  #'1',  # arb. units
        # TODO: FIXME: smooth nacs should be Hartree/Bohr?
        'smooth_nacs': '1',  # arb. units
        'dipoles': traj["dipoles"].attrs["unit"]
        if "dipoles" in traj
        else "1",  #  '1', # arb. units
    }

    if "velocities" in traj:
        metadata['_property_unit_dict']["velocities"] = traj["velocities"].attrs["unit"]

    metadata['atomrefs'] = {}

    metadata['n_singlets'] = traj.attrs["num_singlets"]  # 3  # S0, S1, and S2
    metadata['n_doublets'] = traj.attrs["num_doublets"]
    metadata['n_triplets'] = traj.attrs["num_triplets"]  # 0  # no triplets

    # TODO: FIXME: Not sure if we extract this from the data?
    metadata['phasecorrected'] = (
        False  # phase-properties (NACs, dipoles) are not phase corrected
    )

    metadata['states'] = " ".join(
        traj.state_names.values
    )  # 'S S S'  # three singlet states

    return metadata


@needs(data_vars=set(["energy", "atNames", "atNums", "atXYZ"]))
def write_ase_db(
    traj: Trajectory,
    db_path: str,
    db_format: Literal['schnet', 'spainn'] | None,
    keys_to_write: Collection | None = None,
    preprocess: bool = True,
):
    """Function to write a Dataset into a ASE db in either SchNet or SPaiNN format.

    Args:
        traj (Trajectory): The Dataset to be written to an ASE db style database
        db_path (str): Path to write the database to
        db_format (Literal["schnet", "spainn";] | None): Format of the target database. Used to control order of dimensions in data arrays. Can be either "schnet" or "spainn".
        keys_to_write (Collection | None, optional): Optional parameter to restrict which data variables to . Defaults to None.
        preprocess (bool, optional): _description_. Defaults to True.

    Raises:
        ValueError: If neither `frame` nor `time` dimension is present on the dataset.
        ValueError: If the `db_format` is neither `schnet`, `spainn` nor None

    Notes:
        See `https://spainn-md.readthedocs.io/en/latest/userguide/data_pipeline.html#generate-a-spainn-database` for details on SPaiNN format.
    """
    leading_dim_name: Literal['frame', 'time']

    if "frame" in traj.dims:
        leading_dim_name = 'frame'
    elif "time" in traj.dims:
        leading_dim_name = 'time'
    else:
        raise ValueError(
            "Neither `frame` nor `time` dimension present in dataset. No leading dimension differentiating between frames could be identified."
        )

    if preprocess:
        traj = _prepare_for_write_schnetpack(traj, leading_dim_name)

    statedims = ['state', 'statecomb', 'state_or_statecomb']
    if db_format == 'schnet':
        order = ['frame', *statedims, 'atom', 'direction']
        traj = traj.transpose(*order, missing_dims='ignore')
    elif db_format == 'spainn':
        traj['energy'] = traj['energy'].expand_dims('tmp', axis=1)
        order = ['frame', 'tmp', 'atom', *statedims, 'direction']
        traj = traj.transpose(*order, missing_dims='ignore')
    elif db_format is None:
        # leave the axis orders as they are
        pass
    else:
        raise ValueError(
            f"'db_format' should be one of 'schnet', 'spainn' or None, not '{db_format}'"
        )

    # TODO: FIXME: Do we really want to tabula rasa existing databases?
    if os.path.exists(db_path):
        os.remove(db_path)

    # Restrict, which data variables are written.
    data_var_keys = set([str(x) for x in traj.data_vars.keys()])
    if not keys_to_write:
        keys_to_write = data_var_keys
    else:
        keys_to_write = data_var_keys.intersection(keys_to_write)
    keys_to_write = keys_to_write.difference(['atNames'])

    with connect(db_path, type='db') as db:
        # FIXME: Metadata is only required for SPaiNN, but it seems to me like there is no harm in applying it to SchNarc as well.
        db.metadata = _collect_metadata(traj)
        for i, frame in traj.groupby(leading_dim_name):
            # Remove leading dimension
            frame = frame.squeeze(leading_dim_name)

            # Set a few key parameters from our input parsing functions
            kv_pairs = {}
            if "delta_t" in traj.variables:
                kv_pairs["delta_t"] = traj["delta_t"].values
            else:
                kv_pairs["delta_t"] = traj.attrs["delta_t"]

            if "input_format" in traj.variables:
                kv_pairs["input_format"] = traj["input_format"].values
            else:
                kv_pairs["input_format"] = traj.attrs["input_format"]

            if "input_type" in traj.variables:
                kv_pairs["input_type"] = traj["input_type"].values
            else:
                kv_pairs["input_type"] = traj.attrs["input_type"]

            if "input_format_version" in traj.variables:
                kv_pairs["input_format_version"] = traj["input_format_version"].values
            else:
                kv_pairs["input_format_version"] = traj.attrs["input_format_version"]

            # Actually output the entry
            db.write(
                Atoms(
                    symbols=frame['atNames'].data,
                    positions=frame['atXYZ'],
                    numbers=frame['atNums'],
                    velocities=frame["velocities"] if "velocities" in frame else None,
                ),
                key_value_pairs=kv_pairs,
                data={k: frame[k].data for k in keys_to_write},
            )
