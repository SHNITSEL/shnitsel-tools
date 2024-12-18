import matplotlib.pyplot as plt
import textalloc as ta


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
        # ax.text(scdata['time'][-1], scdata['mean'][-1], sc, c=c, va=vas[sc], ha='right')
        ta.text(
            ax, scdata['time'][-1], scdata['mean'][-1], sc, c=c, va=vas[sc], ha='right'
        )
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


def plot_timeplots(pops, delta_E, fosc_time, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic([['pop'], ['de'], ['ft']], layout='constrained')
        fig.set_size_inches(8.27 / 3, 11.69 / 2)

    plot_pops(pops, axs['pop'])
    plot_time_interstate_error(delta_E, axs['de'])
    if fosc_time is not None:
        plot_time_interstate_error(fosc_time, axs['ft'])

    for axn in ['de', 'pop']:
        axs[axn].sharex(axs['ft'])
        axs[axn].tick_params(axis='x', labelbottom=False)
    axs['ft'].set_xlabel(r'$t$ / fs')  # TODO
    axs['ft'].minorticks_on()

    return axs