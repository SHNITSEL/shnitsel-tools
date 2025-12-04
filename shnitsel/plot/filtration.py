import matplotlib.pyplot as plt

from shnitsel.core.filtration2 import (
    cum_max_quantiles,
    true_upto,
    cum_mask_from_dataset,
)


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
        quantiles['intercept'] = true_upto(quantiles < filtranda['thresholds'], 'time')

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
    return axs


def validity_populations(ds, intersections=True):
    mask = cum_mask_from_dataset(ds)
    if 'thresholds' in mask.coords:
        mask = mask.drop('thresholds')
    if 'good_throughout' in mask.coords:
        mask = mask.drop('good_throughout')
    mask = (
        mask.to_dataset('criterion')
        .assign({'total_population': mask.coords['is_frame']})
        .to_dataarray('criterion')
    )
    counts = mask.sum('trajid')
    means = counts.mean('time')
    if intersections:
        counts = mask.sortby(means, ascending=False).cumprod('criterion').sum('trajid')
    else:
        counts = counts.sortby(means, ascending=False)
    fig, axs = plt.subplots(2, 1)
    fig.set_size_inches(6, 8)
    for criterion in counts.coords['criterion'].data:
        data = counts.sel(criterion=criterion)
        axs[0].plot(data.coords['time'], data, label=criterion)
        axs[1].plot(data.coords['time'], data * data.coords['time'], label=criterion)
    if intersections:
        order = counts.coords['criterion'].data
        labels = [order[0]] + ['AND ' + x for x in order[1:]]
        axs[0].legend(labels)
    else:
        axs[0].legend()
    axs[0].set_ylabel('# trajectories')
    axs[1].set_ylabel('# frames if transected now')
    axs[1].set_xlabel('time / fs')
    return axs