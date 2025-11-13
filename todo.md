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



