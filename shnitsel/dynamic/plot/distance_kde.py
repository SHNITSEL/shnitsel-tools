import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from shnitsel.dynamic import postprocess as P
from .. import pca_biplot as pb
from ..pca_biplot import figax


def fit_kdes(noodle, distances):
    masks = (distances < 3, distances > 5)
    return tuple(stats.gaussian_kde(noodle.sel(frame=m).T) for m in masks)

def eval_kde(kernel, xx, yy):
    Z = kernel.evaluate(np.c_[xx.ravel(), yy.ravel()].T)
    return Z.reshape(xx.shape) / Z.max()

def fit_and_eval_kdes_dist(frames, distances, fineness=500, extension=0.1):
    noodle, hops = P.pca_and_hops(frames)
    c, t = fit_kdes(noodle, distances)
    noodle = noodle.transpose('frame', 'PC')  # required order for the following 3 lines
    means = noodle.mean(dim='frame').values
    mins = noodle.min(dim='frame').values
    mins -= (means - mins) * extension
    maxs = noodle.max(dim='frame').values
    maxs += (maxs - means) * extension
    ls = np.linspace(mins, maxs, num=fineness).T
    xx, yy = np.meshgrid(ls[0], ls[1])
    return xx, yy, eval_kde(c, xx, yy), eval_kde(t, xx, yy)

def plot_kdes(xx, yy, Zcis, Ztrans, levels=None, fig=None, ax=None):
    fig, ax = figax(fig=fig, ax=ax)
    for Z, c in zip([Zcis, Ztrans], ['purple', 'green']):
        ax.contourf(xx, yy, Z, levels=levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=levels, colors=c, linewidths=0.5)

def biplot_distance_time(frames, at1, at2, levels=None):
    """
    Generates a biplot that visualizes PCA projections and kernel density estimates (KDE) 
    of distances between two specified atoms over time.

    Parameters:
    -----------
    frames : xarray.Dataset
        A dataset containing trajectory frames with atomic coordinates.
    at1 : int
        Index of the first atom for distance calculation.
    at2 : int
        Index of the second atom for distance calculation.
    levels : list of float, optional
        Contour levels for the KDE plot. Defaults to [0.08, 1].

    Returns:
    --------
    kde_data : tuple
        The computed KDE data for the atom-atom distance distribution.

    Notes:
    ------
    - Computes Euclidean distances between `at1` and `at2` across all frames.
    - Performs PCA on trajectory data and visualizes clustering of structural changes.
    - Uses kernel density estimation (KDE) to analyze the distance distributions.
    - Produces a figure with PCA projection, cluster analysis, and KDE plots.
    """
  
    if levels is None:
        levels = [0.08, 1]

    distances = P.distance(frames['atXYZ'], at1, at2)

    # prepare layout
    fig, oaxs = plt.subplots(1, 2, layout='constrained', width_ratios=[3, 2])
    fig.set_size_inches(8.27, 11.69 / 3)  # a third of a page, spanning both columns
    gs = oaxs[0].get_subplotspec().get_gridspec()
    for ax in oaxs:
        ax.remove()
    pcasf = fig.add_subfigure(gs[0])
    pcaax = pcasf.subplots(1, 1)
    structsf = fig.add_subfigure(gs[1])
    structaxs = structsf.subplot_mosaic('ab\ncd')

    # prepare data
    kde_data = fit_and_eval_kdes_dist(frames, distances, fineness=100)
    d = pb.pick_clusters(frames, nbins=4)
    
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    
    mol = pb.show_atom_numbers(frames['atXYZ'].isel(frame=0))

    pb.plot_noodleplot(
        *P.pca_and_hops(frames),
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    pb.plot_clusters3(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=mol,
        labels=list('abcd'),
    )
    plot_kdes(*kde_data, levels=levels, ax=pcaax)

    return kde_data

def biplot_distance(frames, at1, at2, levels=None):
    
    if levels is None:
        levels = [0.08, 1]

    distances = P.distance(frames['atXYZ'], at1, at2)

    # prepare layout
    fig, oaxs = plt.subplots(1, 2, layout='constrained', width_ratios=[3, 2])
    fig.set_size_inches(8.27, 11.69 / 3)  # a third of a page, spanning both columns
    gs = oaxs[0].get_subplotspec().get_gridspec()
    for ax in oaxs:
        ax.remove()
    pcasf = fig.add_subfigure(gs[0])
    pcaax = pcasf.subplots(1, 1)
    structsf = fig.add_subfigure(gs[1])
    structaxs = structsf.subplot_mosaic('ab\ncd')

    # prepare data
    kde_data = fit_and_eval_kdes_dist(frames, distances, fineness=100)
    d = pb.pick_clusters(frames, nbins=4)
    
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    
    mol = pb.show_atom_numbers(frames['atXYZ'].isel(frame=0))

    pb.plot_noodleplot(
        *P.pca_and_hops(frames),
        c=distances,
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    pb.plot_clusters3(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=mol,
        labels=list('abcd'),
    )
    
    plot_kdes(*kde_data, levels=levels, ax=pcaax)

    return kde_data
