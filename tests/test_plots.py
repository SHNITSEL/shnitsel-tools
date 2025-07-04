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

FIXDIR = 'tests/fixtures'


@pytest.fixture
def ensembles():
    names = ['butene']
    return {
        name: sh.open_frames(os.path.join(FIXDIR, name, 'data.nc')) for name in names
    }


@pytest.fixture
def spectra3d(ensembles):
    return {
        name: frames.sh.get_inter_state().sh.assign_fosc().sh.spectra_all_times()
        for name, frames in ensembles.items()
    }


#################
# plot.spectra3d:


@image_comparison(['ski_plots'])
def test_ski_plots(spectra3d):
    for name, spectral in spectra3d.items():
        name
        # os.path.join(FIXDIR, name, 'ski_plots.png')
        # with tempfile.NamedTemporaryFile() as f:
        sh.plot.ski_plots(spectral)  # .savefig(f.name)


def test_pcm_plots(): ...


###########
# plot.kde:
def test_biplot_kde(): ...


def test_plot_kdes(): ...


def test_plot_cdf_for_kde(): ...


##############################
# Functions from "pca_biplot":


def test_plot_noodleplot(): ...


def test_plot_noodlelplot_lines(): ...  # once implemented!


def test_plot_loadings(): ...


## Following two together
@pytest.fixture
def highlight_pairs(): ...  # careful -- this uses rdkit, not mpl. What's the return type? Annotate!


def test_mpl_imshow_png(highlight_pairs): ...  # maybe in combination with the above


def test_plot_clusters(): ...  # can we find better names for these? Maybe they're all special cases of a more general function?


def test_plot_clusters2(): ...


def test_plot_clusters3(): ...


def test_plot_bin_edges(): ...


############################
# Functions from "plotting":


def test_pca_line_plot(): ...  # can we generalize this and use the result to finish implementing plot_noodleplot_lines()?


def test_pca_scatter_plot(): ...  # this is unimplemented, and if implemented would be identical to plot_noodleplot, I expect.


def test_timeplot(): ...  # Legacy timeplot function using seaborn via conversion to pandas


def test_timeplot_interstate(): ...  # Legacy timeplot function which does something similar to postprocess.get_inter_state() before plotting


###########################################
# Functions from the "datasheet" hierarchy:

# Skip plot/colormaps.py.
# TODO Skip plot/common.py?
# Skip plot/hist.py?


## plot/__init__.py:
def test_plot_datasheet(): ...  # Warning: this will take long to run -- make optional?


## plot/per_state_hist.py
def test_plot_per_state_histograms(): ...


## plot/dip_trans_hist.py
def test_single_hist(): ...


def test_plot_dip_trans_histograms(): ...


def test_plot_spectra(): ...


def test_plot_separated_spectra_and_hists(): ...  # Monster function! Break up?


## plot/nacs_hist.py
def test_plot_nacs_histograms(): ...


## plot/structure.py


# TODO Why is show_atXYZ deprecated? What has replaced it? The composition of xyz_to_mol() and mol_to_png()?
def test_plot_structure(): ...


## plot/time.py
def test_plot_time_interstate_error(): ...  # TODO 3 statecombs hard-coded for label positioning! Bad!


def test_plot_pops(): ...


def test_plot_timeplots(): ...