import os
import tempfile

import pytest
from matplotlib.testing.decorators import image_comparison

import shnitsel as sh
import shnitsel.xarray

# In this file, we aim to directly test the output of all plotting functions,
# by comparing their output for a test dataset to a pre-made reference plot.
# This does nothing to guarantee the correctness of the reference, but it
# does make it obvious when the graphics are altered by changes to code,
# and when newly-introduced bugs prevent plotting from completing.

# For now we have made the decision not to test the plot-targeting calculation
# backend directly, as this should be subject to thorough change, at least
# before the initial release.

# Framework for now: matplotlib.testing
# Later: matplotcheck (additional dev dependency)

FIXDIR = 'tutorials/test_data/fixtures'



class TestPlotFunctionality:
    """Class to test all plotting functionality included in Shnitsel-Tools"""


    @pytest.fixture
    def ensembles(self):
        names = ['butene']
        return {
            name: sh.open_frames(os.path.join(FIXDIR, name, 'data.nc')) for name in names
        }


    @pytest.fixture
    def spectra3d(self, ensembles):
        return {
            name: frames.sh.get_inter_state().sh.assign_fosc().sh.spectra_all_times()
            for name, frames in ensembles.items()
        }

    #################
    # plot.spectra3d:

    @image_comparison(['ski_plots'])
    def test_ski_plots(self, spectra3d):
        for name, spectral in spectra3d.items():
            name
            # os.path.join(FIXDIR, name, 'ski_plots.png')
            # with tempfile.NamedTemporaryFile() as f:
            sh.plot.ski_plots(spectral)  # .savefig(f.name)

    def test_biplot(self):
        # load trajectory data of A01
        a01 = sh.open_frames('tutorials/test_data/shnitsel/A01_ethene_dynamic.nc')
        # create PCA plot over all trajectories with visualization of the
        # four most important PCA-axis on the molecular structure
        # C=C bond color highlighgting via KDE in PCA
        sh.plot.biplot_kde(
            frames=a01, at1=0, at2=1, geo_filter=[[0, 3], [5, 20]], levels=10
        )
        # C-H bond color highlighting via KDE in PCA
        sh.plot.biplot_kde(
            frames=a01, at1=0, at2=2, geo_filter=[[0, 3], [5, 20]], levels=10
        )

    def test_ski_plots_accessor_conversion(self):
        # load data
        spectra_data = (
            sh.open_frames(path='tutorials/test_data/shnitsel/A01_ethene_dynamic.nc')
            .sh.get_inter_state()
            .sh.spectra_all_times()
        )
        # plot spectra at different simulation times in one plot with a dahsed line that tracks the maximum
        sh.plot.ski_plots(spectra_data)

    def test_pcm_plots(self):
        # TODO
        ...

    ###########
    # plot.kde:
    def test_biplot_kde(self):
        # TODO
        ...

    def test_plot_kdes(self):
        # TODO
        ...

    def test_plot_cdf_for_kde(self):
        # TODO
        ...

    ##############################
    # Functions from "pca_biplot":

    def test_plot_noodleplot(self):
        # TODO
        ...

    def test_plot_noodlelplot_lines(self):  # once implemented!
        # TODO
        ...

    def test_plot_loadings(self):
        # TODO
        ...

    ## Following two together
    @pytest.fixture
    def highlight_pairs(self):
        # careful -- this uses rdkit, not mpl. What's the return type? Annotate!
        # TODO
        ...

    def test_mpl_imshow_png(self,highlight_pairs):
        # maybe in combination with the above
        # TODO
        ...

    def test_plot_clusters(self):
        # can we find better names for these? Maybe they're all special cases of a more general function?
        # TODO
        ...

    def test_plot_clusters2(self):
        # TODO
        ...

    def test_plot_clusters3(self):
        # TODO
        ...

    def test_plot_bin_edges(self):
        # TODO
        ...

    ############################
    # Functions from "plotting":

    def test_pca_line_plot(self):
        # can we generalize this and use the result to finish implementing plot_noodleplot_lines()?

        # TODO
        ...

    def test_pca_scatter_plot(self):
        # this is unimplemented, and if implemented would be identical to plot_noodleplot, I expect.
        # TODO
        ...

    def test_timeplot(self):
        # Legacy timeplot function using seaborn via conversion to pandas

        # TODO
        ...

    def test_timeplot_interstate(self):
        # Legacy timeplot function which does something similar to postprocess.get_inter_state() before plotting

        # TODO
        ...

    ###########################################
    # Functions from the "datasheet" hierarchy:

    # Skip plot/colormaps.py.
    # TODO Skip plot/common.py?
    # Skip plot/hist.py?

    ## plot/__init__.py:
    def test_plot_datasheet(self):
        # Warning: this will take long to run -- make optional?

        # TODO
        ...

    ## plot/per_state_hist.py
    def test_plot_per_state_histograms(self):
        # TODO
        ...

    ## plot/dip_trans_hist.py
    def test_single_hist(self):
        # TODO
        ...

    def test_plot_dip_trans_histograms(self):
        # TODO
        ...

    def test_plot_spectra(self):
        # TODO
        ...

    def test_plot_separated_spectra_and_hists(self):
        # Monster function! Break up?
        # TODO
        ...

    ## plot/nacs_hist.py
    def test_plot_nacs_histograms(self):
        # TODO
        ...

    ## plot/structure.py

    # TODO Why is show_atXYZ deprecated? What has replaced it? The composition of xyz_to_mol() and mol_to_png()?
    def test_plot_structure(self):
        # TODO
        ...

    ## plot/time.py
    def test_plot_time_interstate_error(self):
        # TODO 3 statecombs hard-coded for label positioning! Bad!
        ...

    def test_plot_pops(self):
        # TODO
        ...

    def test_plot_timeplots(self):
        # TODO
        ...
