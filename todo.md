# Missing files

We are currently suffering from missing input files for tests:
- Needs to be replaced: tutorials/test_data/ase/tobias_cis_new.db
- Needs to be replaced: tutorials/test_data/ase/old_CH2NH2.db
- TODO: Add tests for conversion
- Should parsed input-settings be present in the resulting trajectory?
- Check whether NACS are actually available for different formats and what units they have: unit_dimensions.nacs: "1", # TODO: FIXME: NACS in molcas: 1/Bohr, SHARC liest 1/Bohr
- Check reader routines: Several data for energy, forces, nacs, e_kin, phases is full of NANs or the same identical value.
- NACS in SHARC written to output.dat (settings as well?)
- TODO: Mask array to only consider relevant data. Would allow us to limit atom, timestep, etc. Should probably be generated at loading time and not really stored.
- TODO: Make it possible to merge a trajectory into an already multi-trajectory. May be necessary for some users or for us extending datasets in the future.
- TODO: State name assignment is broken. Only state prefix, no state number. Multiple places where states are assigned, should only happen once.
- TODO: FIXME: NewtonX has no state type info, breaking the state names in conversion file. NewtonX has float activestate


