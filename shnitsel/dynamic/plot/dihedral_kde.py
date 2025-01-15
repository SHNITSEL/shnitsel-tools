import numpy as np

from shnitsel.dynamic import postprocess as P
from ..pca_biplot import figax


def fit_kdes(noodle, dihedrals):
    from scipy import stats

    masks = (dihedrals < 80, dihedrals > 110)
    return tuple(stats.gaussian_kde(noodle.sel(frame=m).T) for m in masks)


def eval_kde(kernel, xx, yy):
    Z = kernel.evaluate(np.c_[xx.ravel(), yy.ravel()].T)
    return Z.reshape(xx.shape) / Z.max()


def fit_and_eval_kdes(frames, fineness=500, extension=0.1):
    noodle, hops = P.pca_and_hops(frames)
    dihedrals = P.dihedral(frames['atXYZ'], 0, 1, 2, 3) * 180 / np.pi
    c, t = fit_kdes(noodle, dihedrals)
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
    # if levels is not None:
    #     levels = np.array(levels)

    # cmap = mpl.colormaps['cividis_r']
    # cnorm = mpl.colors.Normalize(0, 1)
    # cscale = mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap)
    # fig.colorbar(cscale, ax=ax)
    for Z, c in zip([Zcis, Ztrans], ['purple', 'green']):
        ax.contourf(xx, yy, Z, levels=levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=levels, colors=c, linewidths=0.5)


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