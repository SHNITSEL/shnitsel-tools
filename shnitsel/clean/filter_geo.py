from dataclasses import dataclass
from typing import Literal, Sequence, TypeVar
from copy import copy

import numpy as np
import xarray as xr
from rdkit.Chem import Mol

from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.filtering.structure_selection import SMARTSstring, StructureSelection
from shnitsel.geo.geocalc import get_distances
from shnitsel.bridges import construct_default_mol
from shnitsel.clean.common import dispatch_filter
from shnitsel.units.conversion import convert_length
from shnitsel.clean.dispatch_plots import dispatch_plots
from shnitsel.units.definitions import length

TrajectoryOrFrames = TypeVar(
    "TrajectoryOrFrames", bound=Trajectory | Frames | xr.Dataset
)


@dataclass
class GeometryFiltrationThresholds:
    """Helper class to keep an extensible set of threshold values for various geometric
    filtration properties.
    The dictionary `match_thresholds` maps SMARTS to a threshold for the length of the bonds described by these SMARTS.
    Additionally, this class has a `length_unit` field which specifies the unit that all other constraints are given in.
    By default, the unit is `angstrom`.

    TODO: FIXME: Should probably provide properties that automatically write to the dict instead.
    """

    all_bonds_threshold: float = 3.0
    all_bonds_smarts: str = "[*]~[*]"
    all_h_to_C_or_N_bonds_threshold: float = 2.0
    all_h_to_C_or_N_bonds_SMARTS: str = "[#6,#7][H]"
    length_unit: str = length.Angstrom

    # Mappings of arbitrary SMARTs to threshold values.
    # Each SMARTs should ideally only cover one bond.
    match_thresholds: dict[str, float] | None = None

    def __init__(self, settings: dict[str, float] | None = None):
        if settings is None:
            self.match_thresholds = dict()
        elif isinstance(settings, dict):
            self.match_thresholds = settings
        elif isinstance(settings, GeometryFiltrationThresholds):
            # initializing from instance of self
            self.match_thresholds = settings.get_full_match_dict()
            self.length_unit = settings.length_unit
        else:
            raise ValueError()

    def get_full_match_dict(self) -> dict[str, float]:
        """Get the full dictionary of SMARTs strings and associated bond length thresholds.

        Used to incorporate settings that are available as explicit fields of this class into the match_thresholds dict.

        Returns:
            dict[str, float]: The combination of settings that can be set via the fields in this class and the `match_thresholds` dict.
        """
        res_dict = dict(self.match_thresholds) if self.match_thresholds else {}
        if self.all_bonds_smarts not in res_dict:
            res_dict[self.all_bonds_smarts] = self.all_bonds_threshold
        if self.all_h_to_C_or_N_bonds_SMARTS not in res_dict:
            res_dict[self.all_h_to_C_or_N_bonds_SMARTS] = (
                self.all_h_to_C_or_N_bonds_threshold
            )

        return res_dict

    def to_dataarray(
        self, selected_SMARTS: Sequence[str] | None = None
    ) -> xr.DataArray:
        """Helper function to convert this dataclass object into an xarray.DataArray to be
        assigned to the coordinate of a Filtration dataset.

        Args:
            selected_criteria (Sequence[str] | None, optional):
                The sequence of criteria keys to be a result of the conversion. Defaults to None.
                Must be a subset of the keys of fields on this object if not set to None.
                If set to None, all but the unit field in this object will be used to set the
                list of criterion keys.


        Returns:
            xr.DataArray:
                A DataArray with a 'criterion' coordinate with either all available or all
                selected properties and the threshold values in the array cells.

        """
        match_dict = self.get_full_match_dict()

        if selected_SMARTS is None:
            selected_SMARTS = list(match_dict.keys())
        threshold_values = [
            match_dict.get(smarts, np.nan) for smarts in selected_SMARTS
        ]

        # Build DataArray from thresholds, criterion names and length unit
        res = xr.DataArray(
            list(threshold_values),
            coords={"criterion": selected_SMARTS},
            attrs={"units": self.length_unit},
        )
        return res.astype(float)


def calculate_bond_length_filtranda(
    frames: xr.Dataset | DataSeries,
    geometry_thresholds: dict[SMARTSstring, float]
    | GeometryFiltrationThresholds
    | None = None,
    mol: Mol | None = None,
) -> xr.DataArray:
    """Derive bond length filtration targets from an xr.Dataset

    Parameters
    ----------
    frames
        A xr.Dataset with an ``atXYZ`` variable
    geometry_thresholds, optional
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
                for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
                for a given search, the longest bond-length will be considered for each frame
            - The unit for the maximum length is provided in the member variable `length_unit` which defaults to `angstrom`.
            - If not provided will be initialized with thresholds for H-(C/N) bonds and one for all bonds.
    mol, optional
        Optional `Mol` object to perform SMARTs matching on.

        TODO: FIXME: We should be able to provide a StructureSelection instead.

    Returns
    -------
        An xr.DataArray of filtration targets stacked along the ``criterion`` dimension;
        one criterion per ``search_dict`` entry.
    """
    if isinstance(frames, xr.Dataset):
        frames = wrap_dataset(frames, DataSeries)

    # Assign default threshold rules.
    if not isinstance(geometry_thresholds, GeometryFiltrationThresholds):
        geometry_thresholds = GeometryFiltrationThresholds(geometry_thresholds)

    thresholds_array = geometry_thresholds.to_dataarray()

    converted_coords = convert_length(
        frames.positions, to=geometry_thresholds.length_unit
    )

    if mol is None:
        mol = construct_default_mol(frames.dataset)
    base_selection = StructureSelection.init_from_mol(mol, ["bonds"])

    criteria_ordered = thresholds_array.coords["criterion"].values
    criteria_results = []

    for smarts in criteria_ordered:
        # Find all bonds conforming to this smarts
        smarts_selection = base_selection.select_bonds(smarts)
        # Calculate distances for these bonds
        smart_specific_distances = get_distances(
            converted_coords, structure_selection=smarts_selection
        )

        # Find maximum across all descriptors/bonds in each frame for this smarts
        max_distances = smart_specific_distances.max("descriptor", keep_attrs=True)
        # Add criterion dimension and append to results
        criteria_results.append(max_distances.expand_dims("criterion"))

    return xr.concat(criteria_results, dim="criterion").assign_coords(
        {"thresholds": thresholds_array}
    )


# TODO: FIXME: This should operate on single trajectories.
def filter_by_length(
    frames_or_trajectory: TrajectoryOrFrames,
    filter_method: Literal["truncate", "omit", "annotate"] | float = "truncate",
    *,
    geometry_thresholds: dict[SMARTSstring, float]
    | GeometryFiltrationThresholds
    | None = None,
    mol: Mol | None = None,
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: Literal["independent", "intersections", False] = False,
) -> TrajectoryOrFrames | None:
    """Filter trajectories according to bond length

    Parameters
    ----------
    frames_or_trajectory: Trajectory | Frames | xr.Dataset
        A Trajectory or Frames Dataset with an ``atXYZ`` variable (NB. this function takes an xr.Dataset as
        opposed to an xr.DataArray for consistency with :py:func:`shnitsel.clean.filter_by_energy`)
    filter_method: Literal["truncate", "omit", "annotate"] | float, optional
        Specifies the manner in which to remove data;

            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
                (:py:func:`shnitsel.clean.truncate`)
            - if 'annotate', merely annotate the data;
            - if a `float` number, interpret this number as a time, and cut all trajectories off at this time,
                discarding those which violate criteria before reaching the given limit,
                (:py:func:`shnitsel.clean.transect`)
        see :py:func:`shnitsel.clean.dispatch_filter`.
    geometry_thresholds: GeometryFiltrationThresholds, optional
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
                for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
                for a given search, the longest bond-length will be considered for each frame
            - The unit for the maximum length is provided in the member variable `length_unit` which defaults to `angstrom`.
            - If not provided will be initialized with thresholds for H-(C/N) bonds and one for all bonds.
    mol: Mol, optional
        An rdkit mol object, if not provided it will be generated from the XYZ coordinates in the data
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations: Literal["independent", "intersections", False], optional
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot

    Returns
    -------
        The filtered Dataset or None if the filter method results in the trajectory being rejected.

    Notes
    -----
    The resulting object has a ``filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten. An existing 'criterion' dimension will be dropped from
    the `frames_or_trajectory` parameter along with all variables and coordinates tied to it.
    """

    analysis_data: Trajectory | Frames = wrap_dataset(
        frames_or_trajectory, Trajectory | Frames
    )

    filtranda = calculate_bond_length_filtranda(
        analysis_data, geometry_thresholds=geometry_thresholds, mol=mol
    )
    frames_dataset = type(analysis_data)(
        analysis_data.dataset.drop_dims(["criterion"], errors="ignore").assign(
            filtranda=filtranda
        )
    )

    dispatch_plots(filtranda, plot_thresholds, plot_populations)

    filter_res = dispatch_filter(frames_dataset, filter_method)

    if isinstance(frames_or_trajectory, xr.Dataset):
        return filter_res.dataset
    else:
        return filter_res
