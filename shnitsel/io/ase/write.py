import json
import os
from typing import Any, Collection, Iterable, Literal

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
        dipoles = traj['dipoles'].data
        dip_attributes = traj['dipoles'].attrs
        # Special case: We already found a dipoles variable
        return traj
    elif 'dip_perm' in traj and 'dip_trans' in traj:
        dip_perm = (
            traj['dip_perm'].transpose(leading_dim_name, 'state', 'direction').data
        )
        dip_trans = (
            traj['dip_trans'].transpose(leading_dim_name, 'statecomb', 'direction').data
        )
        dipoles = np.concat((dip_perm, dip_trans), axis=1)
        dip_attributes = traj['dip_perm'].attrs

        del traj['dip_perm'], traj['dip_trans']
    elif 'dip_perm' in traj:
        dipoles = traj['dip_perm'].data
        dip_attributes = traj['dip_perm'].attrs
        del traj['dip_perm']
    elif 'dip_trans' in traj:
        dipoles = traj['dip_trans']
        dip_attributes = traj['dip_trans'].attrs
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


def _ndarray_to_json_ser(value):
    return {"__ndarray": {"entries": value.tolist(), "dtype": value.dtype.descr}}


def _collect_metadata(traj: Trajectory, keys_to_write: Iterable[str]) -> dict[str, Any]:
    """Helper function to generate the SPaiNN Metadata dict from a Trajectory struct.

    Extracts info from attributes and variables to set up the dict.

    Args:
        traj (Trajectory): The Dataset to extract the metadata from.

    Returns:
        dict[str, Any]: The resulting metadata dictionary.
    """
    # Define metadata information (dictionary)
    metadata: dict[str, Any] = {}
    shnitsel_meta = {}

    if "trajectory_input_path" in traj.attrs:
        metadata['info'] = traj.attrs["trajectory_input_path"]

    if "est_level" in traj.attrs:
        metadata['ReferenceMethod'] = traj.attrs["est_level"]
        # (
        #     'SA3-CASSCF(2,2)'  # state-average CASSCF with 2 electrons in 2 orbitals
        # )

    metadata['_distance_unit'] = (
        traj["atXYZ"].attrs["unit"] if "unit" in traj["atXYZ"].attrs else "Bohr"
    )

    metadata['_property_unit_dict'] = {
        k: (traj[k].attrs["unit"] if "unit" in traj[k].attrs else "1")
        for k in keys_to_write
        if k in traj and "unitdim" in traj[k].attrs
        # 'energy': traj["energy"].attrs["unit"]
        # if "energy" in traj
        # else "1",  # 'Hartree',
        # 'forces': traj["forces"].attrs["unit"]
        # if "forces" in traj
        # else "1",  # 'Hartree/Bohr',
        # 'nacs': traj["nacs"].attrs["unit"]
        # if "nacs" in traj
        # else "1",  #'1',  # arb. units
        # # TODO: FIXME: smooth nacs should be Hartree/Bohr?
        # 'smooth_nacs': '1',  # arb. units
        # 'dipoles': traj["dipoles"].attrs["unit"]
        # if "dipoles" in traj
        # else "1",  #  '1', # arb. units
    }

    if "dipoles" in traj:
        metadata['_property_unit_dict']["dipoles"] = (
            traj["dipoles"].attrs["unit"]
            if "dipoles" in traj and "unit" in traj["dipoles"].attrs
            else "1"
        )

    # if "velocities" in traj:
    #     metadata['_property_unit_dict']["velocities"] = traj["velocities"].attrs["unit"]

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

    # Very specific Shnitsel stuff:
    shnitsel_meta["misc_attrs"] = {
        k: v for k, v in traj.attrs.items() if not str(k).startswith("__")
    }

    shnitsel_meta["var_meta"] = {
        varname: {
            "attrs": {
                v_k: v
                for v_k, v in traj[varname].attrs.items()
                if not str(v_k).startswith("__")
            },
            "dims": [str(d) for d in traj[varname].dims],
        }
        for varname in traj.variables.keys()
    }
    shnitsel_meta["coords"] = {
        coordname: {
            "values": _ndarray_to_json_ser(traj[coordname].values),
            "dims": [str(d) for d in traj[coordname].dims],
        }
        for coordname in traj.coords.keys()
        if coordname not in traj.indexes
        # or traj.indexes[coordname].name != coordname
        or len(traj.indexes[coordname])
        <= 1  # Do not store variables from multi-indices here
    }
    shnitsel_meta["dims"] = {
        dimname: {"length": traj.sizes[dimname]} for dimname in traj.sizes.keys()
    }

    midx_names = []
    for name, index in traj.indexes.items():
        if index.name == name and len(index.names) > 1:
            if len(midx_names) == 0:
                shnitsel_meta["__multi_indices"] = {}
            midx_names.append(name)
            midx_levels = list(index.names)

            shnitsel_meta["__multi_indices"][f'_MultiIndex_levels_for_{name}'] = {
                "level_names": midx_levels,
                "index_tuples": index.values.tolist(),
            }

    shnitsel_meta['_MultiIndex_levels_from_attrs'] = 1

    metadata["__shnitsel_meta"] = shnitsel_meta

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

    statedims = ['state', 'statecomb', 'full_statecomb', 'state_or_statecomb']
    if db_format == 'schnet':
        order = [leading_dim_name, *statedims, 'atom', 'direction']
        traj = traj.transpose(*order, missing_dims='ignore')
    elif db_format == 'spainn':
        traj['energy'] = traj['energy'].expand_dims('tmp', axis=1)
        order = [leading_dim_name, 'tmp', 'atom', *statedims, 'direction']
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
    keys_to_write = keys_to_write.difference(['atNames', 'velocities', 'atXYZ'])

    with connect(db_path, type='db') as db:
        # FIXME: Metadata is only required for SPaiNN, but it seems to me like there is no harm in applying it to SchNarc as well.
        meta_dict = _collect_metadata(traj, keys_to_write)
        meta_dict['n_steps'] = traj.sizes[leading_dim_name]
        if db_format is not None:
            meta_dict["__shnitsel_meta"]["db_format"] = db_format

        meta_dict["__shnitsel_meta"]['shnitsel_leading_dim'] = leading_dim_name
        meta_dict["__shnitsel_meta"] = json.dumps(meta_dict["__shnitsel_meta"])

        db.metadata = meta_dict

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

            if "time" in frame:
                print(frame["time"])
                info = {"time": frame["time"][0]}
            else:
                info = {}

            if "trajid" in frame:
                print(frame["trajid"])
                info["trajid"] = frame["trajid"][0]

            # Actually output the entry
            db.write(
                Atoms(
                    symbols=frame['atNames'].values,
                    positions=frame['atXYZ'].values,
                    # numbers=frame['atNums'],
                    velocities=frame["velocities"] if "velocities" in frame else None,
                    # info={"frame_attrs": info_attrs},
                    info=info,
                ),
                key_value_pairs=kv_pairs,
                data={k: frame[k].data for k in keys_to_write},
            )
