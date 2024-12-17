import matplotlib.pyplot as plt


def plot_time_interstate_error(data, dcol_inter, ylabel, ax):
    vas = {
        '$S_2 - S_0$': 'bottom',
        '$S_2 - S_1$': 'bottom',
        '$S_1 - S_0$': 'top',
    }
    for sc, scdata in data.groupby('statecomb'):
        c = dcol_inter[sc]
        scdata = scdata.squeeze('statecomb')
        ax.fill_between('time', 'upper', 'lower', data=scdata, color=c, alpha=0.3)
        ax.plot('time', 'mean', data=scdata, c=c, lw=0.5)
        ax.text(scdata['time'][-1], scdata['mean'][-1], sc, c=c, va=vas[sc], ha='right')
    ax.set_ylabel(ylabel)
    return ax


def plot_pops(pops, dcol_state, ax):
    for state, sdata in pops.groupby('state'):
        c = dcol_state[state]
        ax.plot(sdata['time'], sdata, c=c, lw=0.5)
        ax.text(sdata['time'][-1], sdata[-1], r"$S_%d$" % state, c=c)
    ax.set_ylabel('Population')
    return ax


def plot_timeplots(pops, delta_E, fosc_time, dcol_state, dcol_inter, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic([['pop'], ['de'], ['ft']], layout='constrained')
        fig.set_size_inches(8.27 / 3, 11.69 / 2)

    plot_pops(pops, dcol_state, axs['pop'])
    plot_time_interstate_error(delta_E, dcol_inter, r'$\Delta E$ / eV', axs['de'])
    if fosc_time is not None:
        plot_time_interstate_error(
            fosc_time, dcol_inter, r'$f_\mathrm{osc}$', axs['ft']
        )

    for axn in ['de', 'pop']:
        axs[axn].sharex(axs['ft'])
        axs[axn].tick_params(axis='x', labelbottom=False)
    axs['ft'].set_xlabel(r'$t$ / fs')
    axs['ft'].minorticks_on()

    return axs