import datetime
from os import PathLike
import numpy as np
from timeit import default_timer as timer
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.figure import Figure

from .datasheet_page import DatasheetPage
from ...data.shnitsel_db.db_function_decorator import concat_subtree
from ...data.shnitsel_db.grouping_methods import group_subtree_by_metadata
from ...data.shnitsel_db_format import ShnitselDB
from ...data.shnitsel_db_helpers import (
    aggregate_xr_over_levels,
    get_trajectories_with_path,
)
from ...data.trajectory_format import Trajectory
from ...io import read

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

_Datasheet_default_page_key = "root"


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

    name: str | None = None
    datasheet_pages: dict[str, DatasheetPage] = {}
    data_source: ShnitselDB | Trajectory

    def __init__(
        self,
        data: Trajectory | ShnitselDB | str | PathLike | Self,
        *,
        name: str | None = None,
        spectra_times: list[int | float] | np.ndarray | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        """Constructor of a datasheet instance.
        If multiple trajectories are provided as a ShnitselDB, a multi-page figure will be generated
        and one page per automatically grouped set of Trajectories will be plotted.



        Args:
            data (Trajectory | ShnitselDB | str | PathLike | Self): Trajectory data as either an individual (possibly concatenated)
                Trajectory object or as a collection of Trajectory objects contained in a ShnitselDB instance.
                Alternatively, a path can be provided from which the data can be loaded via the shnitsel.io.read() function.
                As a last option, another Datasheet instance can be provided and this new instance will be a copy of the other Datasheet.
            name (str, optional): The name of this Datasheet. Will be used as a title for output files if set.
            spectra_times (list[int  |  float] | np.ndarray | None, optional): _description_. Defaults to None.
            col_state (list | None, optional): _description_. Defaults to None.
            col_inter (list | None, optional): _description_. Defaults to None.

        Raises:
            TypeError: If the provided (or read) data is not of Trajectory or ShnitselDB format.

        """
        base_data: Trajectory | ShnitselDB
        self.name = name

        if isinstance(data, Datasheet):
            self._copy_data(old=data)
        else:
            if isinstance(data, str) or isinstance(data, PathLike):
                base_data = read(data, concat_method='db')  # type: ignore # Should be Trajectory or Database
            elif isinstance(data, ShnitselDB) or isinstance(data, Trajectory):
                base_data = data
            else:
                raise TypeError(
                    f"The provided data is neither a Datasheet, a path to Trajectory data or a Trajectory or ShnitselDB object. Was {type(data)}"
                )

            if isinstance(base_data, ShnitselDB):
                self.data_source = base_data
                # TODO: FIXME: Still need to deal with the appropriate grouping of ShnitselDB entries.

                grouped_data = group_subtree_by_metadata(base_data)
                assert (
                    grouped_data is not None and isinstance(grouped_data, ShnitselDB)
                ), "Grouping of the provided ShnitselDB did not yield any result. Please make sure your database is well formed and contains data."

                tree_res_concat = aggregate_xr_over_levels(
                    grouped_data,
                    lambda x: concat_subtree(x, True),
                    "group",
                )
                assert (
                    tree_res_concat is not None
                ), "Aggregation of ShnitselDB yielded None. Please provide a database with data."

                datasheet_groups: list[tuple[str, Trajectory]] = (
                    get_trajectories_with_path(tree_res_concat)
                )

                for name, traj in datasheet_groups:
                    self.datasheet_pages[name] = DatasheetPage(
                        traj,
                        spectra_times=spectra_times,
                        col_inter=col_inter,
                        col_state=col_state,
                    )
                    self.datasheet_pages[name].name = name
            elif isinstance(base_data, Trajectory):
                self.data_source = base_data
                self.datasheet_pages[_Datasheet_default_page_key] = DatasheetPage(
                    self.data_source,
                    spectra_times=spectra_times,
                    col_inter=col_inter,
                    col_state=col_state,
                )
                pass
            else:
                raise TypeError(
                    f"The provided (or read) data is neither Trajectory data nor a Trajectory or ShnitselDB object. Was {type(base_data)}"
                )

    def _copy_data(self, old: Self):
        """Create a copy of an existing Datasheet instance.

        Args:
            old (Self): The old instance to copy
        """
        for key, page in old.datasheet_pages.items():
            self.datasheet_pages[key] = DatasheetPage(page)

        self.data_source = old.data_source
        self.name = old.name

    # @cached_property
    # def axs(self):

    def calc_all(self):
        """Method to precalculate all relevant properties on all (sub-)DatasheetPages"""
        for k, page in self.datasheet_pages.items():
            page.calc_all()

    def plot(
        self,
        include_per_state_hist: bool = False,
        borders: bool = False,
        consistent_lettering: bool = True,
        single_key: str | None = None,
        path: str | PathLike | None = None,
        **kwargs,
    ) -> dict[str, Figure] | Figure:
        """Function to plot datasheets for all trajectory groups/datasets in this Datasheet instance.

        Will output the multi-page figure to a file at `path` if provided.
        Always returns an array of all generated figures to process further.

        Args:
            include_per_state_hist (bool, optional): Flag to include per-state histograms in the plot. Defaults to False.
            borders (bool, optional): A flag whether to draw borders around plots. Defaults to False.
            consistent_lettering (bool, optional): Flag to decide, whether same plots should always have the same letters. Defaults to True.
            single_key (str, optional): Key to a single entry in this set to plot. Keys are specified as paths in the ShnitselDB structure.
            path (str | PathLike | None, optional): Optional path to write a (multi-page) pdf of the resulting datasheets to. Defaults to None.
            **kwargs: Can provide keyword arguments to be used in the pdf metadata dictionary. Among others: 'title', 'author', 'subject', 'keywords'.

        Returns:
            dict[str, Figure]: Map of the keys of the individual datasets to the resulting figure containing all of the Datasheet plots. If no key is available e.g. because a single trajectory was provided, the default key will be "root".
            Figure: If a single_key is specified, will only return that single figure.
        """
        if single_key is None:
            relevant_keys = list(self.datasheet_pages.keys())
        else:
            if single_key not in self.datasheet_pages:
                raise KeyError(
                    f"Provided key {single_key} not found in datasheet pages. Available keys are: {list(self.datasheet_pages.keys())}."
                )
            relevant_keys = [single_key]

        page_figures = {}

        for key in relevant_keys:
            page = self.datasheet_pages[key]
            page_fig = page.plot(
                include_per_state_hist=include_per_state_hist,
                borders=borders,
                consistent_lettering=consistent_lettering,
            )
            page_figures[key] = page_fig

        if path is not None:
            with PdfPages(path) as pdf:
                for key, page_fig in page_figures:
                    pdf.attach_note(f"Plot of: {key}")
                    pdf.savefig(
                        page_fig
                    )  # or you can pass a Figure object to pdf.savefig
            d = pdf.infodict()
            d['Title'] = (
                kwargs['title']
                if 'title' in kwargs
                else (
                    self.name if self.name is not None else 'Shnitsel-Tools Datasheet'
                )
            )
            d['Author'] = kwargs['author'] if 'author' in kwargs else 'Shnitsel-Tools'
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

        if single_key is None:
            return page_figures
        else:
            return page_figures[single_key]

    def _test_subfigures(
        self, include_per_state_hist: bool = False, borders: bool = False
    ):
        """Internal function to test whether subfigure plotting works as intended

        Args:
            include_per_state_hist (bool, optional): Flag to include per-state histograms. Defaults to False.
            borders (bool, optional): Whether the figures should have borders. Defaults to False.
        """
        for key, page in self.datasheet_pages.items():
            page._test_subfigures(
                include_per_state_hist=include_per_state_hist, borders=borders
            )
