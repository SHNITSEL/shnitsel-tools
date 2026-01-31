import json
import logging
import os
from typing import Any, Collection, Iterable, Literal

from ase import Atoms
from ase.db import connect
from ase.db.core import Database
import numpy as np
import xarray as xr

from shnitsel._contracts import needs
from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.tree.node import TreeNode


def _prepare_for_write_schnetpack(
    traj: DataSeries, leading_dim_name: Literal['frame', 'time']
) -> xr.Dataset:
    """Helper function to perform some preprocessing on the dataset before writing to a SchnetPack compatible database.

    Combines the dipole variables into one entry.

    Parameters
    ----------
    traj : DataSeries
        The Dataset to transform into a SchnetPack conforming format.
    leading_dim_name : Literal['frame', 'time']
        The name of the leading dimension identifying different frames within the dataset. Depending on the setup, this should be 'frame' or 'time'.

    Returns
    -------
    xr.Dataset
        The transformed dataset
    """
    # Recombine permanent and transition dipoles, as schnetpack expects
    dipoles: np.ndarray | xr.DataArray | None = None
    working_ds = traj.dataset.copy(deep=False)
    dip_attributes = {}
    if 'dipoles' in traj:
        dipoles = working_ds['dipoles'].data
        dip_attributes = working_ds['dipoles'].attrs
        # Special case: We already found a dipoles variable
        return working_ds
    elif 'dip_perm' in traj and 'dip_trans' in traj:
        dip_perm = (
            working_ds['dip_perm']
            .transpose(leading_dim_name, 'state', 'direction')
            .data
        )
        dip_trans = (
            working_ds['dip_trans']
            .transpose(leading_dim_name, 'statecomb', 'direction')
            .data
        )
        dipoles = np.concat((dip_perm, dip_trans), axis=1)
        dip_attributes = working_ds['dip_perm'].attrs

        del working_ds['dip_perm'], working_ds['dip_trans']
    elif 'dip_perm' in traj:
        dipoles = traj['dip_perm'].data
        dip_attributes = traj['dip_perm'].attrs
        del working_ds['dip_perm']
    elif 'dip_trans' in traj:
        dipoles = traj['dip_trans']
        dip_attributes = traj['dip_trans'].attrs
        del working_ds['dip_trans']

    if dipoles is not None:
        # Change some attributes before assigning
        dip_attributes["long_name"] = "Combined dipole entry (dip_perm and dip_trans)"
        dip_attributes["description"] = (
            "Combined dipole moment containing both permanent and transitional dipole information (if available)"
        )

        working_ds['dipoles'] = (
            [leading_dim_name, 'state_or_statecomb', 'direction'],
            dipoles,
            dip_attributes,
        )

    return working_ds


def _ndarray_to_json_ser(value):
    return {"__ndarray": {"entries": value.tolist(), "dtype": value.dtype}}


def _collect_metadata(
    traj: DataSeries,
    keys_to_write: Collection[str],
) -> dict[str, Any]:
    """Helper function to generate the SPaiNN Metadata dict from a Trajectory struct.

    Extracts info from attributes and variables to set up the dict.
    We expect trees to have been converted to common units before this is invoked.

    Parameters
    ----------
    traj : Trajectory
        The Dataset to extract the metadata from.
    keys_to_write: Collection[str], optional
        The keys of variables to write to the db.

    Returns
    -------
    dict[str, Any]
        The resulting metadata dictionary.
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
            "values": traj[coordname].values.tolist(),
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


# TODO: FIXME: Check the return type and Tree interaction
@needs(data_vars=set(["energy", "atNames", "atNums", "atXYZ"]))
def write_ase_db(
    traj: xr.Dataset | DataSeries | TreeNode[Any, DataSeries | xr.Dataset],
    db_path: str,
    db_format: Literal['schnet', 'spainn'] | None = None,
    keys_to_write: Collection[str] | None = None,
    preprocess: bool = True,
    force: bool = False,
):
    """Function to write a Dataset into a ASE db in either SchNet or SPaiNN format.

    Parameters
    ----------
    traj : xr.Dataset | DataSeries | TreeNode[Any, DataSeries | xr.Dataset]
        The Dataset to be written to an ASE db style database. Can also be in tree format.
        If provided as a tree, the data must be consistent with each other, i.e. all coordinates except for the leading dimension must match.
        Inconsistencies
    db_path : str
        Path to write the database to
    db_format : Literal["schnet", "spainn";] | None, optional
        Format of the target database. Used to control order of dimensions in data arrays. Can be either "schnet" or "spainn".
    keys_to_write : Collection | None, optional
        Optional parameter to restrict which data variables to . Defaults to None.
    preprocess : bool, optional
        Whether to apply preprocessing of the data. Defaults to True.
    force: bool, optional
        A flag to force overwriting of an existing database at the position denoted by `db_path`.

    Raises
    ------
    ValueError
        If neither `frame` nor `time` dimension is present on the dataset.
    ValueError
        If an unsupported data type was provided as an input.
    ValueError
        If the `db_format` is neither `schnet`, `spainn` nor None
    ValueError
        If the data in a provided tree is inconsistent.

    Notes
    -----
    See `https://spainn-md.readthedocs.io/en/latest/userguide/data_pipeline.html#generate-a-spainn-database` for details on SPaiNN format.
    """
    leading_dim_name: Literal['frame', 'time'] | None

    # TODO: FIXME: Do we really want to tabula rasa existing databases?
    if os.path.exists(db_path):
        if force:
            logging.info("Removing database at `{db_path}` before write.")
            os.remove(db_path)
        else:
            msg = "The database at `%s` already exists. To avoid data loss, write will not proceed. If you wish to overwrite the existing databse, please set `force=True` on the call to `write_ase_db()`"
            logging.error(
                msg,
                db_path,
            )
            raise FileExistsError(msg % db_path)
    converted_traj: DataSeries
    if isinstance(traj, TreeNode):
        converted_traj = traj.map_data(
            lambda x: wrap_dataset(x, DataSeries).convert(), dtype=DataSeries
        ).as_stacked  # type: ignore
    else:
        converted_traj = wrap_dataset(traj, DataSeries).convert()

    leading_dim_name: Literal['frame', 'time'] = converted_traj.leading_dimension

    assert leading_dim_name in {'frame', 'time'}, (
        "Neither `frame` nor `time` dimension present in dataset. No leading dimension differentiating between frames could be identified."
    )

    if preprocess:
        working_dataset = _prepare_for_write_schnetpack(
            converted_traj, leading_dim_name
        )
    else:
        working_dataset = converted_traj.dataset

    statedims = ['state', 'statecomb', 'full_statecomb', 'state_or_statecomb']

    if db_format == 'schnet':
        order = [leading_dim_name, ..., *statedims, 'atom', 'direction']
        working_dataset = working_dataset.transpose(*order, missing_dims='ignore')
    elif db_format == 'spainn':
        working_dataset['energy'] = working_dataset['energy'].expand_dims('tmp', axis=1)
        order = [leading_dim_name, ..., 'tmp', 'atom', *statedims, 'direction']
        working_dataset = working_dataset.transpose(*order, missing_dims='ignore')
    elif db_format is None:
        # leave the axis orders as they are
        pass
    else:
        raise ValueError(
            f"'db_format' should be one of 'schnet', 'spainn' or None, not '{db_format}'"
        )

    # Restrict, which data variables are written.
    data_var_keys = set([str(x) for x in working_dataset.data_vars.keys()])
    if not keys_to_write:
        keys_to_write = data_var_keys
    else:
        keys_to_write = data_var_keys.intersection(keys_to_write)
    keys_to_write = keys_to_write.difference(['atNames', 'velocities', 'atXYZ'])

    with connect(db_path, type='db') as db:
        # This performs an implicit type check
        wrapped_input = wrap_dataset(working_dataset, Trajectory | Frames)

        # FIXME: Metadata is only required for SPaiNN, but it seems to me like there is no harm in applying it to SchNarc as well.
        meta_dict = _collect_metadata(wrapped_input, keys_to_write)
        meta_dict['n_steps'] = wrapped_input.sizes[wrapped_input.leading_dim]

        if db_format is not None:
            meta_dict["__shnitsel_meta"]["db_format"] = db_format

        if isinstance(wrapped_input, MultiSeriesDataset):
            meta_dict["__shnitsel_meta"]['multi_set_data'] = True
        else:
            meta_dict["__shnitsel_meta"]['multi_set_data'] = False

        meta_dict["__shnitsel_meta"] = json.dumps(meta_dict["__shnitsel_meta"])

        db.metadata = meta_dict

        _write_trajectory_to_db(
            wrapped_input,
            db,
            # path="/", # Was supposed to help with trees.
            keys_to_write=keys_to_write,
        )


def _write_trajectory_to_db(
    traj: Trajectory | Frames,
    db: Database,
    # path: str,
    keys_to_write: Collection[str],
):
    # Set a few key parameters from our input parsing functions
    kv_pairs = {}
    kv_pairs["charge"] = float(traj.charge)
    kv_pairs["max_ts"] = int(traj.max_timestep)
    kv_pairs["t_max"] = float(traj.t_max)
    kv_pairs["delta_t"] = float(traj.delta_t)
    kv_pairs["input_format"] = traj.input_format
    kv_pairs["input_type"] = traj.input_type
    # NOTE: This guard prefix was introduced, because ASE kept interpreting
    # The 2.0 version string as a string convertible to float and refused to work.
    kv_pairs["input_format_version"] = f"v__{traj.input_format_version}"

    for i, frame in traj.groupby(traj.leading_dim):
        # Remove leading dimension
        frame = frame.squeeze(traj.leading_dim)
        local_kv: dict[str, float | str] = dict(kv_pairs)

        if "time" in frame:
            float_time = float(frame["time"])
            # print(frame["time"], "-->", float_time)
            # info = {"time": float_time}
            local_kv["time"] = float_time
        # else:
        #     info = {}

        if "trajid" in frame:
            int_id = int(frame["trajid"])
        elif "atrajectory" in frame:
            int_id = int(frame["atrajectory"])
        elif "trajectory" in frame:
            int_id = int(frame["trajectory"])
        else:
            int_id = 0
        # info["trajectory"] = int_id  # path + str(int_id)
        local_kv["trajectory"] = int_id

        for coordname in frame.coords:
            if frame.coords[coordname].size == 1:
                if coordname in frame.indexes:
                    index = frame.indexes[coordname]
                    if coordname == index.name and len(index.names):
                        continue

                if coordname not in {
                    'trajid',
                    'trajid_',
                    'trajectory',
                    'atrajectory',
                    'time',
                }:
                    if coordname not in local_kv:
                        coord_data = frame.coords[coordname].item()
                        if isinstance(coord_data, (str, int, float, bool)):
                            local_kv[coordname] = coord_data

        # NOTE: I Contrary to the ASE documentation, the info-dict is not written to DB and therefore not restored
        # upon read.
        # print(info)

        # Actually output the entry
        db.write(
            Atoms(
                symbols=frame['atNames'].values,
                positions=frame['atXYZ'].values,
                # numbers=frame['atNums'],
                velocities=frame["velocities"] if "velocities" in frame else None,
                # info={"frame_attrs": info_attrs},
                # info=info,
            ),
            key_value_pairs=local_kv,
            data={k: frame[k].data for k in keys_to_write},
        )
