# Missing files

We are currently suffering from missing input files for tests:
- TODO: Add tests for conversion
- TODO: Mask array to only consider relevant data. Would allow us to limit atom, timestep, etc. Should probably be generated at loading time and not really stored.
- TODO: Make it possible to merge a trajectory into an already multi-trajectory. May be necessary for some users or for us extending datasets in the future.
- TODO: FIXME: Improve error collection for multi-trajectory parsing.
- TODO: FIXME: Improve errors and warnings texts
- Do NewtonX/PyrAI2md have charge information? SHARC 4.0 has charge info. -> Not that we know of. Probably hidden in QM interface. Add annotation tools for users to change charges on trajectories.
- TODO: Allow setting charge in annotation script? 
- TODO: FIXME: Do not unwrap a single trajectory being read in a multi-trajectory setup. Behavior is quite unexpected.
- TODO: Add benchmark for parallelized loading to compare to sequential loading
- TODO: FIXME: Shape mismatch between DCM in PyRAI2md and default trajectory setup
- TODO: Support reading SHARC netcdf output files

- TODO: Add command line tools for annotation and inspection of trajectories/netcdf files
- TODO: Add tutorial for further CLI tools
- TODO: Add tox for testing shnitsel with different python versions?

- TODO: Implement structure selection methods in structure selection.
- TODO: Rework Tree structure to use lightweight tree class and only use DataTree for load/save. Allows for better type hints and more versatile data storage


python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/sharc/iconds_butene/ -o tutorials/test_data/playground/iconds_butene.nc -c butene -est bravo -basis gulasch -log debug > tutorials/test_data/playground/iconds_butene.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C02_sc-mCy_C4N2H9/init -o tutorials/test_data/playground/C02_sc-mCy_C4N2H9.nc -c C02 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C01_sc-Cy_C3N2H7/init -o tutorials/test_data/playground/C01_sc-Cy_C3N2H7.nc -c C01 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log