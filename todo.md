# Missing files

We are currently suffering from missing input files for tests:
- Needs to be replaced: tutorials/test_data/ase/tobias_cis_new.db
- Needs to be replaced: tutorials/test_data/ase/old_CH2NH2.db
- TODO: Add tests for conversion
- Check whether NACS are actually available for different formats and what units they have: unit_dimensions.nacs: "1", # TODO: FIXME: NACS in molcas: 1/Bohr, SHARC liest 1/Bohr
- NACS in SHARC written to output.dat (settings as well?)
- TODO: Mask array to only consider relevant data. Would allow us to limit atom, timestep, etc. Should probably be generated at loading time and not really stored.
- TODO: Make it possible to merge a trajectory into an already multi-trajectory. May be necessary for some users or for us extending datasets in the future.
- TODO: FIXME: Improve error collection for multi-trajectory parsing.
- TODO: FIXME: Improve errors and warnings texts
- TODO: FIXME: Do NewtonX/PyrAI2md have charge information?
- TODO: FIXME: PyrAI2md add reading of NAC/dcm/SOC



python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/sharc/iconds_butene/ -o tutorials/test_data/playground/iconds_butene.nc -c butene -est bravo -basis gulasch -log debug > tutorials/test_data/playground/iconds_butene.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C02_sc-mCy_C4N2H9/init -o tutorials/test_data/playground/C02_sc-mCy_C4N2H9.nc -c C02 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C01_sc-Cy_C3N2H7/init -o tutorials/test_data/playground/C01_sc-Cy_C3N2H7.nc -c C01 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log