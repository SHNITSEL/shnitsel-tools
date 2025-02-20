from matplotlib.axes import Axes

from .common import figaxs_defaults, centertext


def plot_time_interstate_error(data, ax):
    vas = {
        '$S_2 - S_0$': 'bottom',
        '$S_2 - S_1$': 'bottom',
        '$S_1 - S_0$': 'top',
    }
    for sc, scdata in data.groupby('statecomb'):
        c = scdata['_color'].item()
        scdata = scdata.squeeze('statecomb')
        ax.fill_between('time', 'upper', 'lower', data=scdata, color=c, alpha=0.3)
        ax.plot('time', 'mean', data=scdata, c=c, lw=0.5)
        va = vas.get(sc, 'baseline')
        ax.text(scdata['time'][-1], scdata['mean'][-1], sc, c=c, va=va, ha='right')
    ylabel = data.attrs['tex']
    if u := data.attrs.get('units'):
        ylabel += f" / {u}"
    ax.set_ylabel(ylabel)
    return ax


def plot_pops(pops, ax):
    for state, sdata in pops.groupby('state'):
        c = sdata['_color'].item()
        ax.plot(sdata['time'], sdata, c=c, lw=0.5)
        ax.text(sdata['time'][-1], sdata[-1], r"$S_%d$" % state, c=c)  # TODO
    ax.set_ylabel('Population')
    return ax

@figaxs_defaults(mosaic=[['pop'], ['de'], ['ft']], scale_factors=(1 / 3, 1 / 2))
def plot_timeplots(pops, delta_E, fosc_time, axs=None, fig=None) -> dict[str, Axes]:
    plot_pops(pops, axs['pop'])
    plot_time_interstate_error(delta_E, axs['de'])
    if fosc_time is not None:
        plot_time_interstate_error(fosc_time, axs['ft'])
        lowest_ax = axs['ft']
        higher_axnames = ['de', 'pop']
    else:
        centertext(r"No $\mathbf{\mu}_{ij}$ data", ax=axs['ft'])
        axs['ft'].get_yaxis().set_visible(False)
        axs['ft'].get_xaxis().set_visible(False)
        lowest_ax = axs['de']
        higher_axnames = ['pop']

    lowest_ax.set_xlabel(r'$t$ / fs')  # TODO
    lowest_ax.minorticks_on()

    for axn in higher_axnames:
        axs[axn].sharex(lowest_ax)
        axs[axn].tick_params(axis='x', labelbottom=False)

    return axs