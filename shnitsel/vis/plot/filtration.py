import matplotlib.pyplot as plt

from shnitsel.clean.common import (
    cum_max_quantiles,
    _true_upto,
    _filter_mask_from_criterion_mask,
)


shnitsel_blue = (44 / 255, 62 / 255, 80 / 255)  # '#2c3e50'
shnitsel_yellow = '#C4A000'
shnitsel_magenta = '#7E5273'
text_color = '#fff'
text_backgroundcolor = (0, 0, 0, 0.2)


def check_thresholds(ds_or_da, quantiles=None):
    """Display graphs illustrating
        1. how many trajectories meet each criterion throughout, and
        2. quantiles of cumulative maxima over time for each criterion, indicating at what times a given
           proportion has failed the criterion

    Parameters
    ----------
    ds_or_da
        Data to plot
    quantiles, optional
        quantiles to display and mark on the right-hand graph, by default None

    Returns
    -------
        The matplotlib ``Axes`` object of the plots
    """
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
        quantiles['intercept'] = _true_upto(quantiles < filtranda['thresholds'], 'time')

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
            ax.axhline(threshold, c=shnitsel_yellow)

        for qval, qdata in data.groupby('quantile'):
            qdata = qdata.squeeze(['criterion', 'quantile'])

            ax.fill_between(
                qdata.coords['time'], qdata, fc=(0, 0, 0, 0.2), ec=(0, 0, 0, 0)
            )
            ax.text(qdata['time'][-1], qdata[-1], f"{qval*100} %", va='center', c='k')

            t_icept = qdata['intercept'].item()
            ax.vlines(t_icept, 0, threshold, color=shnitsel_yellow, ls=':')
            ax.text(
                t_icept,
                threshold,
                f"{qval*100} % <{t_icept}",
                ha='right',
                va='center',
                c=text_color,
                backgroundcolor=text_backgroundcolor,
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
            color=shnitsel_blue,
        )
        if 'thresholds' in data.coords:
            threshold = data.coords['thresholds'].item()
            ax.axhline(threshold, c=shnitsel_yellow)
            # ax.text(
            #     0.5,
            #     threshold,
            #     str(threshold),
            #     ha='center',
            #     va='bottom',
            #     c=text_color,
            #     backgroundcolor=text_backgroundcolor,
            # )
            ax.text(
                0.5,
                threshold,
                f"{threshold}\n{data.coords['proportion'].item() * 100:.0f} %",
                ha='center',
                va='center',
                c=text_color,
                backgroundcolor=text_backgroundcolor,
            )

    axs[-1, 0].set_xlabel('cumulative density\nof per-traj maxima')
    axs[-1, 1].set_xlabel('time / fs')
    return axs


def validity_populations(ds_or_da, intersections=True):
    """Display two plots showing
    1. how many trajectories meet criteria (or combinations thereof) up to a given time
    2. how many frames would remain if the ensemble were transected at a given time
    (see :py:func:`shnitsel.clean.transect`)

    Parameters
    ----------
    ds_or_da
        Data to plot
    intersections, optional
        whether to plot intersections of criteria (how many trajectories still meet criterion 1 AND criterion 2)
        or to consider criteria independently

    Returns
    -------
        The matplotlib ``Axes`` object of the plots
    """
    if hasattr(ds_or_da, 'data_vars'):
        filtranda = ds_or_da['filtranda'].copy()
    else:
        filtranda = ds_or_da.copy()
    mask = _filter_mask_from_criterion_mask(filtranda)
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
    axs[1].set_ylabel('# frames if transected at time')
    axs[1].set_xlabel('time / fs')
    return axs
