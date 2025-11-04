# Missing files

We are currently suffering from missing input files for tests:
- Needs to be replaced: tutorials/test_data/ase/tobias_cis_new.db
- Needs to be replaced: tutorials/test_data/ase/old_CH2NH2.db
- Missing NewtonX trajectory at: tutorials/test_data/newtonx/
- Missing Pyrai2MD trajectory at: tutorials/test_data/pyrai2md/
- Possibly wrong file: tutorials/test-data/fixtures/butene/data.nc
- Read version of different input formats
- Harmonize set meta-attributes
- Rework merging of metadata of trajectories
- Apply state names and types set by user.
- Allow input state of SHARC icond files to be read from QM.in and not just QM.log
- Go over the SHARC input format units. 
- In general: We need to check the units and whether conversion is happening. 
- Test whether meta-data is present and consistently named.
- Check whether each format actually sets all of the variables currently defined. 
- - Unset variables not assigned by the routine.