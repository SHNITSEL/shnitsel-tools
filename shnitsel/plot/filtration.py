import matplotlib.pyplot as plt

from shnitsel.core.filtration2 import cum_max_quantiles, true_upto


def check_thresholds(ds_or_da, quantiles=None):
    if hasattr(ds_or_da, 'data_vars'):
        filtranda = ds_or_da['filtranda'].copy()
    else:
        filtranda = ds_or_da.copy()
    quantiles = cum_max_quantiles(ds_or_da, quantiles=quantiles)

    if 'thresholds' in filtranda.coords:
        good_throughout = (
            (filtranda < filtranda['thresholds']).groupby('trajid').all('frame')
        )
        filtranda['proportion'] = (
            good_throughout.sum('trajid') / good_throughout.sizes['trajid']
        )
        quantiles['intercept'] = true_upto(quantiles < filtranda['thresholds'])

    fig, axs = plt.subplots(
        quantiles.sizes['criterion'],
        2,
        sharex='col',
        sharey='row',
        layout='constrained',
        width_ratios=[1, 2],
    )
    fig.set_size_inches(6, 2 * quantiles.sizes['criterion'])
    for (title, data), ax in zip(quantiles.groupby('criterion'), axs[:, 1]):
        if 'thresholds' in data.coords:
            threshold = data.coords['thresholds'].item()
            ax.axhline(threshold, c='r')

        for qval, qdata in data.groupby('quantile'):
            qdata = qdata.squeeze(['criterion', 'quantile'])

            ax.fill_between(
                qdata.coords['time'], qdata, fc=(0, 0, 0, 0.2), ec=(0, 0, 0, 0)
            )
            ax.text(qdata['time'][-1], qdata[-1], f"{qval*100} %", va='center', c='k')

            ##############################
            # x0 = last_time_where(qdata < threshold).item()
            # y0 = qdata.sel(time=x0).item()
            t_icept = qdata['intercept'].item()
            ax.vlines(t_icept, 0, threshold, color='r', ls=':')
            ax.text(
                t_icept,
                threshold,
                f" <{t_icept}",
                ha='center',
                va='bottom',
                c='r',
                rotation='vertical',
                fontsize=6,
            )
            ax.text(
                t_icept,
                threshold,
                f"{qval*100} %",
                ha='right',
                va='top',
                c='r',
                rotation='vertical',
                fontsize=6,
            )
            # print((qdata - threshold))
            ##############################

    for (title, data), ax in zip(filtranda.groupby('criterion'), axs[:, 0]):
        data = data.squeeze('criterion')
        ax.set_ylabel(title)
        ax.hist(
            data.groupby('trajid').max(),
            density=True,
            cumulative=True,
            orientation='horizontal',
            color='b',
        )
        if 'thresholds' in data.coords:
            threshold = data.coords['thresholds'].item()
            ax.axhline(threshold, c='r')
            ax.text(0.5, threshold, str(threshold), ha='center', va='bottom', c='r')
            ax.text(
                0.5,
                threshold,
                f"{data.coords['proportion'].item()*100} %",
                ha='center',
                va='top',
                c='r',
            )

    axs[-1, 0].set_xlabel('cumulative density\nof per-traj maxima')
    axs[-1, 1].set_xlabel('time / fs')