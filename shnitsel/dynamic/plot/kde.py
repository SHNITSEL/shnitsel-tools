from typing import Literal

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from shnitsel.dynamic import postprocess as P
from .. import pca_biplot as pb
from ..pca_biplot import figax


def fit_kdes(noodle, geo_prop, geo_filter):
    kernels = []
    for p1, p2 in geo_filter:
        mask = (p1 < geo_prop) & (geo_prop < p2)
        subset = noodle.sel(frame=mask).T
        if subset.size == 0:
            raise ValueError(f"No points in range {p1} < x < {p2}")
        kernels.append(stats.gaussian_kde(subset))
    return kernels


def eval_kdes(kernels: list, xx, yy):
    xys = np.c_[xx.ravel(), yy.ravel()].T
    Zs = []
    for k in kernels:
        Z = k.evaluate(xys)
        Z = Z.reshape(xx.shape) / Z.max()
        Zs.append(Z)
    return Zs


def get_xx_yy(noodle, fineness=500, extension=0.1):
    means = noodle.mean(dim='frame').values
    mins = noodle.min(dim='frame').values
    mins -= (means - mins) * extension
    maxs = noodle.max(dim='frame').values
    maxs += (maxs - means) * extension
    ls = np.linspace(mins, maxs, num=fineness).T
    xx, yy = np.meshgrid(ls[0], ls[1])
    return xx, yy


def fit_and_eval_kdes(noodle, geo_prop, geo_filter, fineness=500, extension=0.1):
    noodle = noodle.transpose('frame', 'PC')  # required order for the following 3 lines

    xx, yy = get_xx_yy(noodle, fineness=fineness, extension=extension)
    kernels = fit_kdes(noodle, geo_prop, geo_filter)
    return xx, yy, eval_kdes(kernels, xx, yy)


def plot_kdes(xx, yy, Zs, colors=None, levels=None, fig=None, ax=None):
    fig, ax = figax(fig=fig, ax=ax)
    if colors is None:
        if len(Zs) == 2:
            colors = ['purple', 'green']
        else:
            colors = plt.get_cmap('tab10')

    for Z, c in zip(Zs, colors):
        ax.contourf(xx, yy, Z, levels=levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=levels, colors=c, linewidths=0.5)

def biplot_kde(
    frames,
    at1=0,
    at2=1,
    at3=None,
    at4=None,
    geo_filter=None,
    levels=None,
    scatter_color: Literal['time', 'geo'] = 'time',
):
    if scatter_color not in {'time', 'geo'}:
        raise ValueError("`scatter_color` must be 'time' or 'geo'")

    if levels is None:
        levels = [0.08, 1]

    match [at1, at2, at3, at4]: 
        case [at1, at2, None, None]:
            # compute distance between atoms at1 and at2
            geo_prop = P.distance(frames['atXYZ'], at1, at2)
            if not geo_filter:
                geo_filter = [[0,3], [5,100]]
        case [at1, at2, at3, None]:
            # compute angle between vectors at1 - at2 and at2 - at3
            geo_prop = P.angle(frames['atXYZ'], at1, at2, at3, deg=True)
            if not geo_filter:
                geo_filter = [[0,80], [110,180]]
        case [at1, at2, at3, at4]:
            # compute dihedral defined as angle between normals to planes (at1, at2, at3) and (at2, at3, at4)
            geo_prop = P.dihedral(frames['atXYZ'], at1, at2, at3, at4, deg=True)
            if not geo_filter:
                geo_filter = [[0,80], [110,180]]
  
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
    noodle, hops = P.pca_and_hops(frames)
    kde_data = fit_and_eval_kdes(noodle, geo_prop, geo_filter, fineness=100)
    d = pb.pick_clusters(frames, nbins=4)
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    mol = pb.show_atom_numbers(frames['atXYZ'].isel(frame=0))

    if scatter_color == 'time':
        noodleplot_c = None
        noodleplot_cmap = None
    elif scatter_color == 'geo':
        noodleplot_c = geo_prop
        noodleplot_cmap = 'PRGn'
    else:
        assert False

    pb.plot_noodleplot(
        noodle,
        hops,
        c=noodleplot_c,
        cmap=noodleplot_cmap,
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


def plot_cdf_for_kde(z, level, ax=None):
    fig, ax = figax(ax=ax)
    bins, edges, _ = ax.hist(
        z,
        bins=1000,
        range=(0, 1.1 * level),
        cumulative=True,
        density=True,
        histtype='step',
    )
    y = bins[abs(edges - level).argmin()]
    ax.plot([0, level], [y, y], c='r')
    ax.plot([level, level], [0, y], c='r')
    return y
