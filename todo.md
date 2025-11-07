# Missing files

We are currently suffering from missing input files for tests:
- Needs to be replaced: tutorials/test_data/ase/tobias_cis_new.db
- Needs to be replaced: tutorials/test_data/ase/old_CH2NH2.db
- Go over the SHARC input format units. 
- In general: We need to check the units and whether conversion is happening. 
- Should parsed input-settings be present in the resulting trajectory?
- Check whether NACS are actually available for different formats and what units they have: unit_dimensions.nacs: "1", # TODO: FIXME: NACS in molcas: 1/Bohr, SHARC liest 1/Bohr
- Read settings from output.dat in SHARC directories instead of output.log
- Check reader routines: Several data for energy, forces, nacs, e_kin, phases is full of NANs or the same identical value.
- NACS in SHARC written to output.dat (settings as well?)
- TODO: state, atom, direction, satecomb have no values set. Coordinates are dropped in the end because they are not marked as assigned


