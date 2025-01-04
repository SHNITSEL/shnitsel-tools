import numpy as np


# it will be useful to truncate some of the histograms
# this should be noted in the text
# a logarithmic histogram could show outliers... probably not worth it
def trunc_max(data, rel_cutoff=0.01, bins=1000):
    freqs, edges = np.histogram(
        data, bins=bins, range=(np.nanmin(data), np.nanmax(data))
    )
    cutoff = freqs.max() * rel_cutoff
    relevant = edges[:-1][freqs > cutoff]
    return relevant.max()


def truncate(data, rel_cutoff=0.01, bins=1000):
    sup = trunc_max(data, rel_cutoff=rel_cutoff, bins=bins)
    plot_data = data[data <= sup]
    return plot_data


def create_marginals(ax):
    axx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
    axy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
    # no labels next to main plot, no axis at all on other side
    axx.tick_params(axis="x", labelbottom=False)
    axy.tick_params(axis="y", labelleft=False)
    axx.get_yaxis().set_visible(False)
    axy.get_xaxis().set_visible(False)
    return axx, axy


def create_marginals_dict(axs, label):
    ax = axs[label]
    axx, axy = create_marginals(ax)
    axs[f'{label}x'], axs[f'{label}y'] = axx, axy
    return axx, axy