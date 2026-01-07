import os

import matplotlib.pyplot as plt
import numpy as np

import pytest
from pytest import fixture
from matplotlib.testing.decorators import image_comparison

import shnitsel as sh
import shnitsel.xarray
from shnitsel.data.tree import tree_to_frames

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
            ('tutorials/tut_data/traj_I02.nc', -1),
        ]
    )
    def ensembles(self, request):
        path, charge = request.param
        db = read(path)
        res = tree_to_frames(db)
        res['atXYZ'].attrs['charge'] = charge
        res.attrs['charge'] = charge
        return res

    @pytest.fixture
    def spectra3d(self, ensembles):
        return ensembles.st.get_inter_state().st.assign_fosc().st.spectra_all_times()

    #################
    # plot.spectra3d:

    # @image_comparison(['ski_plots'])
    def test_ski_plots(self, spectra3d):
        sh.vis.plot.ski_plots(spectra3d)

    def test_pcm_plots(self, spectra3d):
        sh.vis.plot.spectra3d.pcm_plots(spectra3d)

    ###########
    # plot.kde:
    def test_biplot_kde(self, ensembles):
        sh.vis.plot.kde.biplot_kde(
            frames=ensembles,
            at1=0,
            at2=1,
            geo_filter=[[0.0, 20.0]],
            levels=10,
        )

    @pytest.fixture
    def kde_data(self, ensembles):
        noodle, _ = sh.analyze.pca.pca_and_hops(ensembles, mean=False)
        geo_prop = np.zeros(noodle.sizes['frame'])
        return sh.vis.plot.kde.fit_and_eval_kdes(noodle, geo_prop, [(-1, 1)])

    def test_plot_kdes(self, kde_data):
        sh.vis.plot.kde.plot_kdes(*kde_data)

    def test_plot_cdf_for_kde(self, kde_data):
        xx, yy, Zs = kde_data
        sh.vis.plot.kde.plot_cdf_for_kde(Zs[0].ravel(), 0.1)

    ##############################
    # Functions from "pca_biplot":

    def test_plot_noodleplot(self, ensembles):
        noodle, hops = sh.analyze.pca.pca_and_hops(ensembles, mean=False)
        sh.vis.plot.pca_biplot.plot_noodleplot(noodle, hops)

    @pytest.fixture
    def clusters_loadings_mols(self, ensembles):
        import shnitsel.xarray

        pca = ensembles['atXYZ'].st.pwdists().st.pca('atomcomb')
        loadings = pca.st.loadings
        clusters = sh.vis.plot.pca_biplot.cluster_loadings(loadings)
        mol = sh.bridges.default_mol(ensembles['atXYZ'])
        return clusters, loadings, mol

    def test_plot_loadings(self, clusters_loadings_mols):
        _, loadings, _ = clusters_loadings_mols
        _, ax = plt.subplots(1, 1)
        sh.vis.plot.pca_biplot.plot_loadings(ax, loadings)

    @pytest.fixture
    def highlight_pairs(self, ensembles):
        mol = sh.bridges.default_mol(ensembles['atXYZ'])
        return sh.rd.highlight_pairs(mol, [(0, 1)])

    def test_mpl_imshow_png(self, highlight_pairs):
        _, ax = plt.subplots(1, 1)
        sh.vis.plot.common.mpl_imshow_png(ax, highlight_pairs)

    def test_plot_clusters(self, clusters_loadings_mols):
        clusters, loadings, _ = clusters_loadings_mols
        sh.vis.plot.pca_biplot.plot_clusters(loadings, clusters)

    def test_plot_clusters_insets(self, clusters_loadings_mols):
        clusters, loadings, mol = clusters_loadings_mols
        _, ax = plt.subplots(1, 1)
        sh.vis.plot.pca_biplot.plot_clusters_insets(ax, loadings, clusters, mol)

    def test_plot_clusters_grid(self, clusters_loadings_mols):
        _, ax = plt.subplots(1, 1)
        clusters, loadings, mol = clusters_loadings_mols
        sh.vis.plot.pca_biplot.plot_clusters_grid(loadings, clusters, mol=mol)

    def test_plot_bin_edges(self, ensembles):
        nbins = 5
        data = sh.vis.plot.pca_biplot.pick_clusters(ensembles, nbins=nbins)
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        sh.vis.plot.pca_biplot.plot_bin_edges(
            data['angles'],
            data['radii'],
            data['bins'],
            data['edges'],
            data['picks'],
            ax,
            range(nbins),
        )

    # NB. Functions from the "datasheet" hierarchy tested in separate file