import os

import xarray as xr

from shnitsel.io import read
from shnitsel.io.sharc.trajectory import read_traj
from shnitsel.io.sharc.initial_conditions import dir_of_iconds
from shnitsel.units import standard_shnitsel_units
from shnitsel.units.definitions import unit_dimensions


class TestSHARC:
    """Class to test functionality related to SHARC initial conditions and trajectories"""

    def test_read_dirs_of_iconds(self):
        path = os.path.join("tutorials", "test_data", "sharc", "iconds_butene")
        iconds = dir_of_iconds(path)

        assert isinstance(iconds, xr.Dataset)

        assert "atXYZ" in iconds
        assert "state" in iconds
        assert "energy" in iconds
        assert "forces" in iconds

        assert 'units' in iconds['atXYZ'].attr and iconds['atXYZ'].attr[
            'units'] == standard_shnitsel_units[unit_dimensions.length]
        assert 'units' in iconds['forces'].attr and iconds['forces'].attr[
            'units'] == standard_shnitsel_units[unit_dimensions.force]
        assert 'units' in iconds['energy'].attr and iconds['energy'].attr[
            'energy'] == standard_shnitsel_units[unit_dimensions.energy]

        for var in iconds:
            if 'unitdim' in iconds[var].attrs:
                assert 'units' in iconds[var].attrs
                unit_dim = iconds[var].attrs['unitdim']
                assert iconds[var].attrs['units'] == standard_shnitsel_units[unit_dim]

    def test_sharc_process_and_output_iconds(self):
        # parse icond data
        iconds_butene = dir_of_iconds(
            path="tutorials/test_data/sharc/iconds_butene"
        )
        # save the parsed data to h5netcdf file
        savepath = os.path.join(os.getcwd(), "tutorials",
                                "test_data", "sharc", "iconds_butene.hdf5")
        iconds_butene.reset_index("statecomb").to_netcdf(
            savepath, engine="h5netcdf")

    def test_read_sharc_trajs_direct(self):
        # parse trajectory data from SHARC output files
        traj_frames_butene = read_traj(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002"
        )

        assert isinstance(traj_frames_butene, xr.Dataset)

    def test_read_sharc_glob_match(self):
        # Read a set of trajectories from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/",
            kind='sharc'
        )
        assert isinstance(traj_frames_butene, xr.Dataset)

    def test_read_sharc_wrapper_direct(self):
        # Directly read one single trajectory from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002",
            kind='sharc'
        )

        assert isinstance(traj_frames_butene, xr.Dataset)
