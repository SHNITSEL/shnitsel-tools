# General TODO:

## Primary TODO:

- I/O:
  - [x] Add option to guess the charge with RDKit if not set
  - [ ] Add support for writing complex trees to ASE (currently only stackable trees supported)
- Tree improvement
  - [ ] Tree function wrappers with tree parameters for non-tree parameters in unwrapped function
  - [ ] Renaming/restructuring support of tree structure
  - [ ] Selection functions on the tree (`compound=?`, `group=?`, `<grouping_param>=?`)
    - [x] Added stubs for `.sel()` and `.isel()`
    - [x] Add implementation capturing appropriate dimensions/indexes.
    - [x] Propagate remaining indexers to data entries.
    - [ ] TODO: Make `trajectory` select the desired entries in `Stacked` sets (should select for both `atrajectory` and `trajectory`)
    - [ ] TODO: Make `trajectory` drop the NaN frames in result for `Layered` sets.
  - [ ] Pattern matching in the `getitem()` method, i.e. `db['/I01/**/data']` `db["/I01/*/{1-20}"]`
  - [ ] TODO: FIXME: typing.get_origin() does not work on TreeNode specializations. Maybe just patch the `base` class in `_create_extended_node_class`?
  - [ ] TODO: Add backup visualization probably with networkx to circumvent issue with html repr not working on github. (See TODO: note in TreeNode class after `_repr_html()`)
- (De)serialization:
  - [x] Implementation of Supports(To/From)XrConversion for various (wrapper) types
  - [ ] Tutorial on how to add own types
- Wrapper types:
  - [ ] Add supported Shnitsel-tools functions as direct methods on wrapper types
  - [ ] Improve Visualization/helper text
  - [ ] Add DataArray wrapper
    - [ ] Add specifically an `AtXYZ`/`Positions` Wrapper type.
  - [ ] Add more wrapper types for support function returns
- StructureSelection
  - [x] Options in function signatures (provide directly a SMARTs string, deal with trees, etc.)
    - [x] On-the-go construction added to `geo._get_default_selection()`
    - [x] Patch existing functions to express support for descriptors
  - [x] Add merge/subtract/intersect operations
  - [x] Add 'BLA' and 'pwdist' feature descriptors.
  - [ ] Fix non-redundant coordinates
  - [ ] Structure selection (raise error if empty/warning if empty)
  - [ ] draw: Draw grid of highlighted feature levels
- [ ] Analogs tree, structure selection, warning if no match for compounds
  - [x] Adapted to tree support
  - [ ] Copy if multiple matches
  - [ ] For now: error if multiple matches
- StateSelection
  - [x] Add Support for textual representation of state selection
  - [x] Add merge/subtract/intersect operations
- biplot_kde needs to be fixed to use the descriptors of the PCA in the side plots.
  - [ ] Fix PCA loadings main contribution plot not being the same as the explained PCA.
  - [ ] Fix cluster loadings plot occasionally empty.
- Visualization support
  - [ ] Add generic `plot()` function to various types
  - [ ] Add option for plots from tree hierarchies
- [ ] Add tutorial for further CLI tools
- Dimension reduction:
  - [ ] Refactor PLS
  - [ ] Refactor LDA
- [ ] Clustering support
  - [ ] DBSCAN sounds like a good option to support
- Datasheet
  - [ ] Improve documentation
    - [x] Refactor documentation
    - [ ] Make docstrings more detailed
  - DatasheetPage:
    - [ ] Improve PCA Page
    - [x] Improve default settings for datasheet and pages.
- Hops
  - [ ] Fix some functions not yet moved to tree types or new style of type selection.
    - [x] `assign_hop_time()` done
    - [ ] `focus_hops()`
    - [ ] PLS
    - [ ] LDA
  
- [ ] CLI tutorial
- [ ] Full tutorial retinal
  - [x] Retinal BLA/HOOP/Dihedrals (Put all dihedrals in), 3->2 Dimension reduction
  - [x] Retinal PCA biplot to show we can find the right parameters
  - [x] Length threshold statistics plot
  - [x] Datasheet (note, that generation takes somewhat longer)
  - [ ] add one more cool plot, the hop-time alignment would be perfect

- st accessors
  - [x] Fix union types not being unwrapped in generator
  - [x] FQTN types not being passed as TypeForms/string for types not explicitly imported.
  - [ ] Replace type vars by their bounds or `Any`.
  - [ ] TODO: FIXME: Generic Parameters get lost in generation.
  

## Secondary TODO:
- [ ] Add tests for conversion
- [ ] Mask array to only consider relevant data. Would allow us to limit atom, timestep, etc. Should probably be generated at loading time and not really stored.
- [ ] Make it possible to merge a trajectory into an already multi-trajectory. May be necessary for some users or for us extending datasets in the future.
- Do NewtonX/PyrAI2md have charge information? SHARC 4.0 has charge info. -> Not that we know of. Probably hidden in QM interface.

- [ ] Add command line tools for annotation and inspection of trajectories/netcdf files
- [ ] Add tox for testing shnitsel with different python versions

## Tertiary TODO:

- [ ] Profiling of key functions like `get_bats()`
- [ ] FIXME: Shape mismatch between DCM in PyRAI2md and default trajectory setup
- [ ] Support reading SHARC netcdf output files

### Profiling notes:
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/sharc/iconds_butene/ -o tutorials/test_data/playground/iconds_butene.nc -c butene -est bravo -basis gulasch -log debug > tutorials/test_data/playground/iconds_butene.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C02_sc-mCy_C4N2H9/init -o tutorials/test_data/playground/C02_sc-mCy_C4N2H9.nc -c C02 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log
python -m cProfile ./shnitsel/cli/convert_to_shnitsel.py tutorials/test_data/playground/C01_sc-Cy_C3N2H7/init -o tutorials/test_data/playground/C01_sc-Cy_C3N2H7.nc -c C01 -est bravo -basis gulasch -log debug > tutorials/test_data/playground/C02.profile.log