import os

import matplotlib.pyplot as plt
import numpy as np

import pytest
from pytest import fixture
from matplotlib.testing.decorators import image_comparison

import shnitsel as st
from shnitsel.analyze.pca import pca_and_hops
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.tree.node import TreeNode
from shnitsel.vis.plot.pca_biplot import cluster_loadings, plot_loadings
import shnitsel.xarray
from shnitsel.bridges import default_mol

from shnitsel.io import read

# In this file, we aim to directly test the output of all plotting functions,
# by comparing their output for a test dataset to a pre-made reference plot.
# This does nothing to guarantee the correctness of the reference, but it
# does make it obvious when the graphics are altered by changes to code,
# and when newly-introduced bugs prevent plotting from completing.

# Framework for now: matplotlib.testing
# Later: matplotcheck (additional dev dependency)


class TestPlotFunctionality:
    """Class to test all plotting functionality included in Shnitsel-Tools"""

    @fixture(
        params=[
            ('./tutorials/test_data/shnitsel/traj_I02.nc', 1),
        ]
    )
    def ensembles(self, request) -> Frames:
        path, charge = request.param
        db = read(path)
        assert isinstance(db, TreeNode)
        res = db.set_charge(charge).as_stacked
        return res

    @pytest.fixture
    def spectra3d(self, ensembles):
        return ensembles.st.get_inter_state().st.assign_fosc().st.spectra_all_times()

    #################
    # plot.spectra3d:

    # @image_comparison(['ski_plots'])
    @pytest.mark.xfail
    def test_ski_plots(self, spectra3d):
        from shnitsel.vis.plot import ski_plots

        ski_plots(spectra3d)

    @pytest.mark.xfail
    def test_pcm_plots(self, spectra3d):
        from shnitsel.vis.plot.spectra3d import pcm_plots

        pcm_plots(spectra3d)

    ###########
    # plot.kde:
    @pytest.mark.xfail
    def test_biplot_kde(self, ensembles):
        from shnitsel.vis.plot.kde import biplot_kde

        biplot_kde(
            ensembles,
            0,
            1,
            geo_kde_ranges=[(0.0, 20.0)],
            contour_levels=10,
        )

    @pytest.fixture
    def kde_data(self, ensembles):
        from shnitsel.vis.plot.kde import _fit_and_eval_kdes

        noodle, _ = pca_and_hops(ensembles, center_mean=False)
        geo_prop = np.zeros(noodle.projected_inputs.sizes['frame'])
        return _fit_and_eval_kdes(noodle, geo_prop, [(-1, 1)])
    
    @pytest.mark.xfail
    def test_plot_kdes(self, kde_data):
        from shnitsel.vis.plot.kde import _plot_kdes

        _plot_kdes(*kde_data)

    @pytest.mark.xfail
    def test_plot_cdf_for_kde(self, kde_data):
        from shnitsel.vis.plot.kde import plot_cdf_for_kde

        xx, yy, Zs = kde_data
        plot_cdf_for_kde(Zs[0].ravel(), 0.1)

    ##############################
    # Functions from "pca_biplot":

    def test_plot_noodleplot(self, ensembles):
        from shnitsel.vis.plot.pca_biplot import plot_noodleplot

        noodle, hops = pca_and_hops(ensembles, center_mean=False)
        plot_noodleplot(noodle.projected_inputs, hops)

    @pytest.fixture
    def clusters_loadings_mols(self, ensembles: Frames):
        import shnitsel.xarray
        from shnitsel.analyze.pca import pca, PCAResult

        pca_res: PCAResult = pca(ensembles)  # ['atXYZ'].st.pwdists().st.pca('atomcomb')
        loadings = pca_res.principal_components
        clusters = cluster_loadings(loadings)
        mol = default_mol(ensembles)
        return clusters, loadings, mol
    
    @pytest.mark.xfail
    def test_plot_loadings(self, clusters_loadings_mols):
        _, loadings, _ = clusters_loadings_mols
        _, ax = plt.subplots(1, 1)
        plot_loadings(ax, loadings)

    @pytest.fixture
    def highlight_pairs(self, ensembles):
        from shnitsel.rd import highlight_pairs

        mol = default_mol(ensembles)
        return highlight_pairs(mol, [(0, 1)])

    @pytest.mark.xfail
    def test_mpl_imshow_png(self, highlight_pairs):
        from shnitsel.vis.plot.common import mpl_imshow_png

        _, ax = plt.subplots(1, 1)
        mpl_imshow_png(ax, highlight_pairs)

    @pytest.mark.xfail
    def test_plot_clusters(self, clusters_loadings_mols):
        from shnitsel.vis.plot.pca_biplot import plot_clusters

        clusters, loadings, _ = clusters_loadings_mols
        plot_clusters(loadings, clusters)

    @pytest.mark.xfail
    def test_plot_clusters_insets(self, clusters_loadings_mols):
        clusters, loadings, mol = clusters_loadings_mols
        from shnitsel.vis.plot.pca_biplot import plot_clusters_insets

        _, ax = plt.subplots(1, 1)
        plot_clusters_insets(ax, loadings, clusters, mol)

    @pytest.mark.xfail
    def test_plot_clusters_grid(self, clusters_loadings_mols):
        from shnitsel.vis.plot.pca_biplot import plot_clusters_grid

        _, ax = plt.subplots(1, 1)
        clusters, loadings, mol = clusters_loadings_mols
        plot_clusters_grid(loadings, clusters, mol=mol)

    @pytest.mark.xfail
    def test_plot_bin_edges(self, ensembles):
        from shnitsel.vis.plot.pca_biplot import plot_bin_edges, pick_clusters

        nbins = 5
        data = pick_clusters(ensembles, num_bins=nbins)
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        plot_bin_edges(
            data['angles'],
            data['radii'],
            data['bins'],
            data['edges'],
            data['picks'],
            ax,
            range(nbins),
        )

    # NB. Functions from the "datasheet" hierarchy tested in separate file
