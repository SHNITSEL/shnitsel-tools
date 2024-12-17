import matplotlib.pyplot as plt
import numpy as np

from .hist import truncate


def plot_per_state_histograms(per_state, dcol_state, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic(
            [['energies', 'forces', 'dip_perm']], layout='constrained'
        )
        fig.set_size_inches(8.27, 11.69 / 5)

    for quantity in ['energies', 'forces', 'dip_perm']:
        ax = axs[quantity]

        for state, data in per_state.groupby('state'):
            c = dcol_state[state - 1]
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

    eunits = (
        per_state['energies'].attrs.get('units') or per_state['energies'].attrs['unit']
    )
    axs['energies'].set_xlabel(rf'$E$ / {eunits}')
    axs['forces'].set_xlabel(r'$\mathbf{F}$ / (hartree/bohr)')
    axs['dip_perm'].set_xlabel(r'$\mathbf{\mu}_i$ / (??)')