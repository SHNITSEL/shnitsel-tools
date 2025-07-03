import numpy as np

from .common import figaxs_defaults, centertext, symbols
from .hist import truncate

@figaxs_defaults(mosaic=[['energy', 'forces', 'dip_perm']], scale_factors=(1, 1 / 5))
def plot_per_state_histograms(per_state, axs=None, fig=None):
    for quantity in ['energy', 'forces', 'dip_perm']:
        ax = axs[quantity]
        if quantity not in per_state:
            centertext("No %s data" % symbols.get(quantity, quantity), ax)
            continue

        for state, data in per_state.groupby('state'):
            c = data['_color'].item()
            counts, edges, _ = ax.hist(
                truncate(data[quantity].squeeze(), bins=100),
                color=c,
                alpha=0.2,
                bins=100,
            )
            ax.plot((edges[1:] + edges[:-1]) / 2, counts, c=c, lw=0.5)
            idxmax = np.argmax(counts)
            ax.text(
                edges[[idxmax, idxmax + 1]].mean(),
                counts[idxmax],
                r"$S_%d$" % (state - 1),
                c=c,
            )

        long_name = per_state[quantity].attrs.get('long_name')
        units = per_state[quantity].attrs.get('units')
        axs[quantity].set_xlabel(rf'{long_name} / {units}')

    # for quantity in ['forces', 'dip_perm']:
    #     axs[quantity].
    axs['energy'].set_ylabel('# points')
    return axs