from dataclasses import dataclass
import datetime
import logging
from os import PathLike
from typing import Any, Mapping, Sequence
from matplotlib.axes import Axes
import numpy as np
from timeit import default_timer as timer
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.figure import Figure, SubFigure
from tqdm import tqdm

from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.data_group import DataGroup
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.state_selection import StateSelection, StateSelectionDescriptor
from shnitsel.filtering.structure_selection import (
    StructureSelection,
    StructureSelectionDescriptor,
)

from .datasheet_page import DatasheetPage
from ...data.shnitsel_db_helpers import concat_subtree
from ...data.tree import ShnitselDB
from ...data.shnitsel_db_helpers import (
    aggregate_xr_over_levels,
)
from ...data.dataset_containers import Trajectory, Frames, wrap_dataset
from ...io import read
import xarray as xr

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

_Datasheet_default_page_key = "root"


@dataclass
class Datasheet:
    """Class to generate overview plots for a collection of trajectories.

    Multiple individual plots are possible.
    Available plots include:
    - per_state_histograms: Histograms of energy, forces and transition dipoles per state
    - separated_spectra_and_hists: Histograms of transition dipoles and time plots
    - noodle: Noodle plots of structure over time for each states
    - structure: Plot of the moleculare structure given either all positions or a smiles map
    - nacs_histograms: A histogram of the nacs between states as well as energy and force histograms
    - timeplots: Plot of the active states over time.
    """

    data_source: (
        TreeNode[Any, ShnitselDataset | xr.Dataset] | ShnitselDataset | xr.Dataset
    )
    datasheet_pages: dict[str, DatasheetPage]
    name: str | None = None

    def __init__(
        self,
        data: ShnitselDataset
        | xr.Dataset
        | ShnitselDB[ShnitselDataset | xr.Dataset]
        | str
        | PathLike
        | Self,
        state_selection: StateSelection | StateSelectionDescriptor | None = None,
        structure_selection: StructureSelection
        | StructureSelectionDescriptor
        | None = None,
        *,
        name: str | None = None,
        spectra_times: list[int | float] | np.ndarray | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        """Constructor of a datasheet instance.
        If multiple trajectories are provided as a ShnitselDB, a multi-page figure will be generated
        and one page per automatically grouped set of Trajectories will be plotted.

        Parameters
        ----------
        data : ShnitselDataset | xr.Dataset | TreeNode[Any, ShnitselDataset | xr.Dataset] | str | PathLike | Self)
            Trajectory data as either an individual (possibly concatenated)
            Trajectory object or as a collection of Trajectory objects contained in a ShnitselDB instance.
            Alternatively, a path can be provided from which the data can be loaded via the shnitsel.io.read() function.
            As a last option, another Datasheet instance can be provided and this new instance will be a copy of the other Datasheet.
        state_selection : StateSelection | StateSelectionDescriptor, optional
            Optional parameter to specify a subset of states and state combinations that may be considered for the dataset.
            Will be generated if not provided.
        structure_selection: StructureSelection | StructureSelectionDescriptor, optional
            Optional parameter to limit the PCA plot and analysis to a specific subset of the structure.
            Will be generated if not provided.
        name : str, optional
            The name of this Datasheet.
            Will be used as a title for output files if set.
        spectra_times : list[int  |  float] | np.ndarray | None, optional
            Sequence of times to calculate spectra at. Defaults to None.
        col_state : list | None, optional
            A list of colors to use for the states. Defaults to default shnitsel colors.
        col_inter : list | None, optional
            A list of colors to use for state combinations. Defaults to default shnitsel colors.

        Raises
        ------
        TypeError
            If the provided (or read) data is not of Trajectory or ShnitselDB format.

        """
        base_data: ShnitselDataset | xr.Dataset | TreeNode
        self.name = name
        self.datasheet_pages = {}

        if isinstance(data, Datasheet):
            self._copy_data(old=data)
        else:
            if isinstance(data, str) or isinstance(data, PathLike):
                base_data = read(data, concat_method='db')  # type: ignore # Should be Trajectory or Database
            elif isinstance(data, TreeNode):
                base_data = data.map_data(wrap_dataset)
            elif isinstance(data, ShnitselDataset):
                base_data = data
            else:
                try:
                    base_data = wrap_dataset(data)
                except:
                    raise TypeError(
                        f"The provided data is neither a Datasheet, a path to Trajectory data or a Trajectory or ShnitselDB object. Was {type(data)}"
                    )

            if isinstance(base_data, TreeNode):
                self.data_source = base_data
                # TODO: FIXME: Still need to deal with the appropriate grouping of ShnitselDB entries.

                grouped_data = base_data.group_data_by_metadata()
                assert grouped_data is not None and isinstance(
                    grouped_data, ShnitselDB
                ), (
                    "Grouping of the provided ShnitselDB did not yield any result. Please make sure your database is well formed and contains data."
                )

                tree_res_concat: TreeNode[Any, ShnitselDataset | xr.Dataset] | None = (
                    grouped_data.map_filtered_nodes(
                        lambda n: isinstance(n, DataGroup) and n.is_flat_group,
                        lambda n: n.construct_copy(
                            children={'agg': DataLeaf(data=n.as_stacked)}
                        ),
                        dtype=Frames,
                    )
                )
                assert tree_res_concat is not None, (
                    "Aggregation of ShnitselDB yielded None. Please provide a database with data."
                )

                datasheet_groups: list[tuple[str, Frames | xr.Dataset]] = list(
                    tree_res_concat.map_filtered_nodes(
                        lambda n: n.is_leaf,
                        lambda n: DataLeaf(name="tmp", data=(n.path, n.data)),
                        tuple,
                    ).collect_data()
                )

                for name, traj in datasheet_groups:
                    self.datasheet_pages[name] = DatasheetPage(
                        traj,
                        spectra_times=spectra_times,
                        state_selection=state_selection,
                        structure_selection=structure_selection,
                        col_inter=col_inter,
                        col_state=col_state,
                    )
                    self.datasheet_pages[name].name = name
            elif isinstance(base_data, (Trajectory, Frames, xr.Dataset)):
                self.data_source = base_data
                self.datasheet_pages[_Datasheet_default_page_key] = DatasheetPage(
                    self.data_source,
                    spectra_times=spectra_times,
                    structure_selection=structure_selection,
                    col_inter=col_inter,
                    col_state=col_state,
                )
                pass
            else:
                raise TypeError(
                    f"The provided (or read) data is neither Trajectory data nor a Trajectory/Frames or ShnitselDB object. Was {type(base_data)}"
                )

    def _copy_data(self, old: Self):
        """Create a copy of an existing Datasheet instance.

        Parameters
        ----------
        old : Self
            The old instance to copy
        """
        for key, page in old.datasheet_pages.items():
            self.datasheet_pages[key] = DatasheetPage(page)

        self.data_source = old.data_source
        self.name = old.name

    # @cached_property
    # def axs(self):

    def calc_all(self):
        """Method to precalculate all relevant properties on all (sub-)DatasheetPages"""
        for page in self.datasheet_pages.values():
            page.calc_all()

    @property
    def pages(self) -> Mapping[str, DatasheetPage]:
        """Retrieve the mapping of individual pages contained in this datasheet.

        Returns
        -------
        Mapping[str, DatasheetPage]
            The keys of the pages are the individual paths in a hierarchical structure,
            whereas the values are the `DataSheetPage` instances.
        """
        return self.datasheet_pages

    def plot(
        self,
        include_per_state_hist: bool = False,
        include_coupling_page: bool = True,
        include_pca_page: bool = False,
        include_meta_page: bool = False,
        borders: bool = False,
        consistent_lettering: bool = True,
        single_key: str | None = None,
        path: str | PathLike | None = None,
        **kwargs,
    ) -> dict[str, list[Figure]] | list[Figure]:
        """Function to plot datasheets for all trajectory groups/datasets in this Datasheet instance.

        Will output the multi-page figure to a file at `path` if provided.
        Always returns an array of all generated figures to process further.

        Parameters
        ----------
        include_per_state_hist : bool, optional
            Flag to include per-state histograms in the plot. Defaults to False.
        include_coupling_page : bool, optional
            Flag to create a full page with state-coupling plots. Defaults to False.
        include_pca_page : bool, optional
            Flag to create a PCA analysis page with details on PCA results. Defaults to True.
        include_meta_page : bool, optional
            Flag to add a page with meta-information about the trajectory data. Defaults to False
        borders : bool, optional
            A flag whether to draw borders around plots. Defaults to False.
        consistent_lettering : bool, optional
            Flag to decide, whether same plots should always have the same letters. Defaults to True.
        single_key : str, optional
            Key to a single entry in this set to plot. Keys are specified as paths in the ShnitselDB structure.
        path : str | PathLike | None, optional
            Optional path to write a (multi-page) pdf of the resulting datasheets to. Defaults to None.
        **kwargs
            Can provide keyword arguments to be used in the pdf metadata dictionary.
            Among others: 'title', 'author', 'subject', 'keywords'.

        Returns
        -------
        dict[str, Figure]
            Map of the keys of the individual datasets to the resulting figure containing all of the Datasheet plots.
            If no key is available e.g. because a single trajectory was provided, the default key will be "root".
        Figure
            If a `single_key` is specified, will only return that single figure.
        """
        if single_key is None:
            relevant_keys = list(self.datasheet_pages.keys())
        else:
            if single_key not in self.datasheet_pages:
                raise KeyError(
                    f"Provided key {single_key} not found in datasheet pages. Available keys are: {list(self.datasheet_pages.keys())}."
                )
            relevant_keys = [single_key]

        page_figures: dict[str, list[Figure]] = {}

        for key in relevant_keys:
            page = self.datasheet_pages[key]
            page_fig = page.plot(
                include_per_state_hist=include_per_state_hist,
                include_coupling_page=include_coupling_page,
                include_pca_page=include_pca_page,
                include_meta_page=include_meta_page,
                borders=borders,
                consistent_lettering=consistent_lettering,
            )
            page_figures[key] = page_fig if isinstance(page_fig, list) else [page_fig]

        if path is not None:
            logging.info(
                "Saving datasheet as pdf. Please be patient, this may take some time."
            )
            with PdfPages(path) as pdf:
                for key, page_fig in tqdm(
                    page_figures.items(), unit="page", desc="Written"
                ):
                    pdf.attach_note(f"Plot of: {key}")
                    for fig in page_fig:
                        pdf.savefig(fig, dpi=300)
                d = pdf.infodict()
                d['Title'] = (
                    kwargs['title']
                    if 'title' in kwargs
                    else (
                        self.name
                        if self.name is not None
                        else 'Shnitsel-Tools Datasheet'
                    )
                )
                d['Creator'] = 'Shnitsel-Tools package'
                d['Author'] = (
                    kwargs['author'] if 'author' in kwargs else 'Shnitsel-Tools'
                )
                d['Subject'] = (
                    kwargs['subject']
                    if 'subject' in kwargs
                    else 'Visualization of key statistics'
                )
                d['Keywords'] = (
                    kwargs['keywords']
                    if 'keywords' in kwargs
                    else 'Datasheet shnitsel shnitsel-tools'
                )
                d['CreationDate'] = (
                    datetime.datetime.today()
                )  # datetime.datetime(2009, 11, 13)
                d['ModDate'] = datetime.datetime.today()
                logging.info(f"Writing pdf with {pdf.get_pagecount()} pages")

        if single_key is None:
            return page_figures
        else:
            return page_figures[single_key]

    def _test_subfigures(
        self, include_per_state_hist: bool = False, borders: bool = False
    ):
        """Internal function to test whether subfigure plotting works as intended

        Parameters
        ----------
        include_per_state_hist : bool, optional
            Flag to include per-state histograms. Defaults to False.
        borders : bool, optional
            Whether the figures should have borders. Defaults to False.
        """
        for page in self.datasheet_pages.values():
            page._test_subfigures(
                include_per_state_hist=include_per_state_hist, borders=borders
            )

    def plot_per_state_histograms(
        self,
        shape: tuple[int, int] | None = None,
    ) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_per_state_histograms()` on each page.

        Parameters
        ----------
        shape : tuple[int, int] | None, optional
            The desired grid shape as `(rows,columns)`, by default None

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_per_state_histograms(
                state_selection=page.state_selection, shape=shape
            )
            for page in self.datasheet_pages.values()
        ]

    def plot_timeplots(self) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_timeplots()` on each page.

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_timeplots(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]

    def plot_separated_spectra_and_hists(
        self,
        current_multiplicity: int = 1,
    ) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_per_state_histograms()` on each page.

        Parameters
        ----------
        current_multiplicity : int, default = 1
            The target multiplicity for which the spectra and histograms should be plotted, by default None

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_separated_spectra_and_hists(
                state_selection=page.state_selection,
                current_multiplicity=current_multiplicity,
            )
            for page in self.datasheet_pages.values()
        ]

    def plot_separated_spectra_and_hists_groundstate(
        self,
    ) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_separated_spectra_and_hists_groundstate()` on each page.

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_separated_spectra_and_hists_groundstate(
                state_selection=page.state_selection
            )
            for page in self.datasheet_pages.values()
        ]

    def plot_nacs_histograms(
        self,
    ) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_nacs_histograms()` on each page.

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_nacs_histograms(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]

    def plot_noodle(
        self,
    ) -> Sequence[Axes]:
        """Helper method to get the results of the call to `plot_noodle()` on each page.

        Returns
        -------
        Sequence[Axes]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_noodle(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]

    def plot_energy_bands(
        self,
    ) -> Sequence[dict[str, Axes]]:
        """Helper method to get the results of the call to `plot_energy_bands()` on each page.

                Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_energy_bands(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]

    def plot_structure(
        self,
    ) -> Sequence[Axes]:
        """Helper method to get the results of the call to `plot_structure()` on each page.

        Parameters
        ----------
        shape : tuple[int, int] | None, optional
            The desired grid shape as `(rows,columns)`, by default None

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_structure(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]

    def plot_pca_structure(
        self,
    ) -> Sequence[Axes]:
        """Helper method to get the results of the call to `plot_pca_structure()` on each page.

        Returns
        -------
        Sequence[dict[str, Axes]]
            The resulting plots per datasheet page.
        """
        return [
            page.plot_pca_structure(state_selection=page.state_selection)
            for page in self.datasheet_pages.values()
        ]
