import logging
from itertools import product
from dataclasses import dataclass
import numpy as np
import scipy.stats as st
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import rdkit

from . import(
  postprocess as P,
  xrhelpers as xh,
  pca_biplot)

from .pca_biplot import plot_noodleplot

# %% [markdown]
# # Calculations

# %% [markdown]
# ## Spectra

# %%
def get_spectrum(data, t, sc, cutoff=0.01):
    data = data.sel(time=t, statecomb=sc)
    res = P.broaden_gauss(data.energies, data.fosc, width=0.1, agg_dim='trajid')
    non_negligible = res.where(res > cutoff, drop=True).energy
    if len(non_negligible) == 0:
        return res.sel(energy=non_negligible)
    return res.sel(
      energy=slice(non_negligible.min(), non_negligible.max()))

def calc_spectra(spectral, times=None, cutoff=0.01):
    """Returns a `dict` of DataArrays indexed by `(time, statecomb)` tuples. 
    """
    if times is None:
        times = [0,10, 20, 30]
    return {
      (t, sc): get_spectrum(spectral, t, sc, cutoff=cutoff)
      for t, sc in product(times, spectral.statecomb.values)
    }

# %%

def get_sgroups(spectra):
    ground, excited = {}, {}
    for (t, sc), v in spectra.items():
        if sc == '$S_2 - S_1$':
            excited[t, sc] = v
        else:
            ground[t, sc] = v
            
    sgroups = (ground, excited)
    return sgroups

# %% [markdown]
# ## Calculate all at once
# TODO: Consider making a dataclass with all the different
# Datasets and methods to process and plot them.

# %%
def calc_all(frames):
    per_state = P.get_per_state(frames)
    inter_state = P.get_inter_state(frames)
    pops = P.calc_pops(frames)

    delta_E = P.time_grouped_ci(inter_state['energies'])
    noodle, hops = P.pca_and_hops(frames)

    if 'dip_trans' in frames:
        inter_state = P.assign_fosc(inter_state)
        fosc_time = P.time_grouped_ci(inter_state['fosc'])
        spectra = calc_spectra(inter_state)
        sgroups = get_sgroups(spectra)
    else:
        fosc_time = None
        spectra = None
        sgroups = None
    
    return {
        'per_state': per_state,
        'inter_state': inter_state,
        'pops': pops,
        'sgroups': sgroups,
        'noodle': noodle,
        'hops': hops,
        'delta_E': delta_E,
        'fosc_time': fosc_time,
        'atXYZ': frames['atXYZ'].isel(frame=0)
    }

# %% [markdown]
# # Helper functions for plotting

clmagma = mpl.colormaps["magma_r"](np.linspace(0, 1, 128))
clmagma[:, 2] *= 1/np.max(clmagma[:, 2]) # more blue near zero, so white rather than yellow
magma_rw = mpl.colors.LinearSegmentedColormap.from_list('magma_rw', clmagma)

# %%
# it will be useful to truncate some of the histograms
# this should be noted in the text
# a logarithmic histogram could show outliers... probably not worth it
def trunc_max(
  data,
  rel_cutoff=0.01,
  bins=1000
):
    freqs, edges = np.histogram(data, bins=bins, range=(np.nanmin(data), np.nanmax(data)))
    cutoff = freqs.max() * rel_cutoff
    relevant = edges[:-1][freqs > cutoff]
    return relevant.max()

def truncate(
  data,
  rel_cutoff=0.01,
  bins=1000
):
    sup = trunc_max(data, rel_cutoff=rel_cutoff, bins=bins)
    plot_data = data[data <= sup]
    return plot_data

def scatter_hist(x, y, ax, ax_histx, ax_histy, main_plotter=None, bins=200, color=None):
    if main_plotter is None:
        main_plotter = lambda x, y: ax.hist2d(x, y, bins=bins, cmap='magma_r')
    
    # the scatter plot:
    main_plotter(x, y)

    # now determine nice limits by hand:
    # binwidth = 0.25
    # xymax = max(np.max(x)), np.max(y)
    # lim = (int(xymax/binwidth) + 1) * binwidth

    # bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, color=color, bins=bins)
    ax_histy.hist(y, orientation='horizontal', color=color, bins=bins)

# Prepare for 2D histograms

def create_marginals(ax):
    axx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
    axy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
    # no labels next to main plot, no axis at all on other side
    axx.tick_params(axis="x", labelbottom=False)
    axy.tick_params(axis="y", labelleft=False)
    axx.get_yaxis().set_visible(False)
    axy.get_xaxis().set_visible(False)
    return axx, axy

# def create_marginals(ax):
#     d = make_axes_locatable(ax)
#     axx = d.append_axes("top", size="10%", pad="4%")
#     axy = d.append_axes("right", size="10%", pad="4%")
#     axx.sharex(ax)
#     axy.sharey(ax)
#     axx.tick_params(axis="x", labelbottom=False)
#     axy.tick_params(axis="y", labelleft=False)
#     axx.get_yaxis().set_visible(False)
#     axy.get_xaxis().set_visible(False)
#     return axx, axy

def create_marginals_dict(axs, label):
    ax = axs[label]
    axx, axy = create_marginals(ax)
    axs[f'{label}x'], axs[f'{label}y'] = axx, axy
    return axx, axy


# %%
# Per-state Histograms

def plot_per_state_histograms(per_state, dcol_state, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic(
            [['energies', 'forces', 'dip_perm']],
            layout='constrained'
        )
        fig.set_size_inches(8.27, 11.69/5)

    for quantity in ['energies', 'forces', 'dip_perm']:
        ax = axs[quantity]
        freqs = []
        edges = []
        for state, data in per_state.groupby('state'):
            c = dcol_state[state-1]
            hd = ax.hist(truncate(data[quantity].squeeze(), bins=100), color=c, alpha=0.2, bins=100)
            # DEBUG:
            if state == 1:
                ret = hd
            counts, edges, _ = hd
            ax.plot((edges[1:] + edges[:-1])/2, counts, c=c, lw=0.5)
            idxmax = np.argmax(counts)
            ax.text(edges[[idxmax, idxmax+1]].mean(), counts[idxmax],
                    r"$S_%d$"%(state-1), c=c)

    eunits = per_state['energies'].attrs.get('units') or per_state['energies'].attrs['unit']
    axs['energies'].set_xlabel(rf'$E$ / {eunits}')
    axs['forces'].set_xlabel(r'$\mathbf{F}$ / (hartree/bohr)')
    axs['dip_perm'].set_xlabel(r'$\mathbf{\mu}_i$ / (??)')

    return ret

    # cruft to be removed soon
        #     f, e, _ = ax.hist(data[quantity], color=c, alpha=0.2, bins=100)
        #     freqs = np.hstack([freqs, f])
        #     edges = np.hstack([edges, e[:-1]]) # one more edge than there are bins
        # cutoff = freqs.max()/100
        # relevant = edges[freqs > cutoff]
        # ax.set_xlim(0, relevant.max())

    # axs['forces'].hist(data['forces'], color=c, alpha=0.2, bins=100)
    # axs['dip_perm'].hist(data['dip_perm'], color=c, alpha=0.2, bins=100)

# %%
# Timeplots

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
        ax.text(scdata['time'][-1], scdata['mean'][-1],
                sc, c=c, va=vas[sc], ha='right')
    ax.set_ylabel(ylabel)
    return ax

def plot_pops(pops, dcol_state, ax):
    for state, sdata in pops.groupby('state'):
        c=dcol_state[state]
        ax.plot(sdata['time'], sdata, c=c, lw=0.5)
        ax.text(sdata['time'][-1], sdata[-1],
                r"$S_%d$"%state, c=c)
    ax.set_ylabel('Population')
    return ax

def plot_timeplots(pops, delta_E, fosc_time, dcol_state, dcol_inter, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic(
          [['pop'],
           ['de'],
           ['ft']],
          layout='constrained')
        fig.set_size_inches(8.27/3, 11.69/2)

    plot_pops(pops, dcol_state, axs['pop'])
    plot_time_interstate_error(delta_E, dcol_inter, r'$\Delta E$ / eV', axs['de'])
    if fosc_time is not None:
        plot_time_interstate_error(fosc_time, dcol_inter, r'$f_\mathrm{osc}$', axs['ft'])

    for axn in ['de', 'pop']:
        axs[axn].sharex(axs['ft'])
        axs[axn].tick_params(axis='x', labelbottom=False)
    axs['ft'].set_xlabel(r'$t$ / fs')
    axs['ft'].minorticks_on()

    return  axs

# %%

clmagma = mpl.colormaps["magma_r"](np.linspace(0, 1, 128))
clmagma[:, 2] *= 1/np.max(clmagma[:, 2]) # more blue near zero, so white rather than yellow
magma_rw = mpl.colors.LinearSegmentedColormap.from_list('magma_rw', clmagma)

custom_ylgnr = mpl.colors.LinearSegmentedColormap.from_list(
      'custom', mpl.colormaps['YlGn_r'](np.linspace(0, 0.75, 128)))

def single_hist(data, shi, slo, color, bins=100, ax=None, cmap=None, cnorm=None):
    if ax is None: _, ax = plt.subplots(1,1)
    if cmap is None: cmap = magma_rw

    axx, axy = create_marginals(ax)
    xdata = data['energies'].squeeze()
    ydata = data['dip_trans'].squeeze()
    xmax = trunc_max(xdata)
    ymax = trunc_max(ydata)
    axx.hist(xdata, range=(0, xmax), color=color, bins=bins)
    axy.hist(ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins)
    hist2d_output = ax.hist2d(
      xdata, ydata, range=[(0, xmax), (0, ymax)],
      bins=bins, cmap=cmap, norm=cnorm)
    
    ax.set_ylabel(r"$\|\mathbf{\mu}_{%d,%d}\|_2$" % (shi, slo))
    ax.text(1.05, 1.05, "$S_%d/S_%d$" % (shi, slo), transform=ax.transAxes,
      ha="left", va="bottom",  color=color,
    #   fontweight='bold',
    )

    return hist2d_output

def plot_dip_trans_histograms(inter_state, dcol_inter, axs=None, cnorm=None):
    if axs is None:
        nplots = len(inter_state.coords['statecomb'])
        _, axs = plt.subplots(nplots, 1, layout='constrained')

    # TODO obviate following cludge:
    sclabels = [(int(x[3]), int(x[9])) for x in inter_state.statecomb.values]

    hist2d_outputs = []
    for i, (sc, data) in enumerate(inter_state.groupby('statecomb')):
        # label = f't{i}'
        shi, slo = sclabels[i]
        ax = axs[i]

        color = dcol_inter[sc]
        hist2d_outputs.append(single_hist(data, shi, slo, color=color, ax=ax, cnorm=cnorm))
    return hist2d_outputs

# %%

def plot_spectra(spectra, dcol_inter, ax=None, cmap=None, cnorm=None):
    if ax is None: _, ax = plt.subplots(1,1)
    cmap = plt.get_cmap(cmap) if cmap else custom_ylgnr
    times = [t for (t, sc) in spectra]
    cnorm = cnorm if cnorm else plt.Normalize(min(times), max(times))
    ax.set_ylabel(r'$f_\mathrm{osc}$')
    ax.invert_xaxis()
    # linestyles = {t: ['-', '--', '-.', ':'][i]
    #               for i, t in enumerate(np.unique(list(zip(*spectra.keys()))[0]))}
    for (t, sc), data in spectra.items():
        c=cmap(cnorm(t))
        # ax.fill_between(data['energy'], data, alpha=0.5, color=c)
        ax.plot(data['energy'], data,
                # linestyle=linestyles[t], c=dcol_inter[sc],
                c=c,
                linewidth=0.5)
    

    return ax

def sep_ground_excited_spectra(spectra, excited_transitions=None):
    if excited_transitions is None:
        excited_transitions = {'$S_2 - S_1$'}

    ground, excited = {}, {}

    for (t, sc), v in spectra.items():
        if sc in excited_transitions:
            excited[t, sc] = v
        else:
            ground[t, sc] = v
    
    return ground, excited

# %%

def _hist_min_max(inter_state):
    minmax = []
    debug = []
    for sc, data in inter_state.groupby('statecomb'):
        data = data.squeeze()
        H, xe, ye = np.histogram2d(data.energies.values, data.dip_trans.values, bins=100)
        minmax += [H.min(), H.max()]
        debug.append((H, xe, ye))
    
    return min(minmax), max(minmax), debug

def plot_separated_spectra_and_hists(
  inter_state, sgroups, dcol_inter, axs=None):
    if axs is None:
        mosaic = [['sg'],
                  ['t0'],
                  ['t1'],
                  ['se'],
                  ['t2'],
                  ['cb_spec'],
                  ['cb_hist']]
        fig, axs = plt.subplot_mosaic(
          mosaic,
          layout='constrained',
          height_ratios=([1]*5)+([0.1]*2))
        fig.set_size_inches(8.27/3,11.69*0.8) # portrait A4

    ground, excited = sgroups
    times = [tup[0] for lst in sgroups for tup in lst]
    scnorm = plt.Normalize(min(times), max(times))
    scmap = plt.get_cmap('viridis_r')
    scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)

    hist2d_outputs = []
    # ground-state spectra and histograms
    plot_spectra(ground, dcol_inter, ax=axs['sg'], cnorm=scnorm, cmap=scmap)
    hist2d_outputs += plot_dip_trans_histograms(
      inter_state.isel(statecomb=[0, 1]),
      dcol_inter,
      axs=[axs[k] for k in ['t1', 't0']])
    # excited-state spectra and histograms
    plot_spectra(excited, dcol_inter, ax=axs['se'], cnorm=scnorm, cmap=scmap)
    hist2d_outputs += plot_dip_trans_histograms(
      inter_state.isel(statecomb=[2]),
      dcol_inter,
      axs=[axs['t2']])

    hists = np.array([tup[0] for tup in hist2d_outputs])
    hcnorm = plt.Normalize(hists.min(), hists.max())

    quadmeshes = [tup[3] for tup in hist2d_outputs]
    for quadmesh in quadmeshes:
        quadmesh.set_norm(hcnorm)

    def ev2nm(ev):
        return 4.135667696*2.99792458*100 / np.where(ev!=0, ev, 1)

    lims = [l for ax in axs.values() for l in ax.get_xlim()]
    new_lims = (min(lims), max(lims))
    for lax, ax in axs.items():
        if lax.startswith('cb'):
            continue
        ax.set_xlim(*new_lims)
        ax.invert_xaxis()

    for ax in list(axs.values()):
        ax.tick_params(axis="x", labelbottom=False)
    axs['t2'].tick_params(axis="x", labelbottom=True)

    secax = axs['sg'].secondary_xaxis('top', functions=(ev2nm, ev2nm))
    secax.set_xticks([50, 75, 100, 125, 150, 200, 300, 500, 1000])
    secax.tick_params(axis='x', rotation=45, labelsize='small')
    for l in secax.get_xticklabels():
        l.set_horizontalalignment('left')
        l.set_verticalalignment('bottom')
    secax.set_xlabel(r'$\Delta E$ / nm')

    for lax in ['cb_spec', 'cb_hist']:
        axs[lax].get_yaxis().set_visible(False)

    axs['cb_spec'].figure.colorbar(scscale, cax=axs['cb_spec'], location='bottom')
    axs['cb_spec'].set_xlabel('time / fs')
    hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
    axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')
    axs['cb_hist'].set_xlabel('# data points')

    axs['se'].set_title(r"$\uparrow$ground state" + "\n" + r"$\downarrow$excited state absorption")
    axs['t2'].set_xlabel(r'$\Delta E$ / eV')

    return axs
        
# %%

# Plot 2D histograms of NACS vs delta_E or dip_trans 
def plot_nacs_histograms(inter_state, hop_idxs, col_inter, axs=None):
    if axs is None:
        fig, axs = plt.subplot_mosaic(
          [['energies', 'forces', 'dip_perm'],
           ['fde', 'pca', 'pca'],
           ['t0', 'pca', 'pca'],
           ['t1', '.', 'pop'],
           ['t2', 'ntd', 'de'],
           ['.', 'nde', 'ft']],
          layout='constrained'
        )
        fig.set_size_inches(8.27,11.69) # portrait A4
    else:
        fig = axs['nde'].figure

    nacs_data = inter_state.sel(frame=hop_idxs)
    axs['nde'].set_ylabel(r'$\Delta E$ / eV')
    axs['nde'].minorticks_on()
    axs['nde'].tick_params(axis="x", labelbottom=False)
    if 'dip_trans' in inter_state:
        axs['ntd'].set_ylabel(r"$\|\mathbf{\mu}_{i,j}\|_2$")
        axs['ntd'].set_xlabel(r"$\|\mathrm{NAC}_{i,j}\|_2$")
        axs['ntd'].minorticks_on()

    def plot(label, yname):
        ax = axs[label]
        axx, axy = create_marginals(ax)
        bins = 100
        # for sc, data in inter_state.groupby('statecomb'):
        for i, (sc, data) in enumerate(nacs_data.groupby('statecomb')):
            if i!=0: continue
            ydata = data[yname].squeeze()
            xdata = data['nacs'].squeeze()
            xmax = trunc_max(xdata)
            ymax = trunc_max(ydata)
            # ymax = trunc_max(ydata) # if you truncate nacs, there's nothing left
            color = col_inter[i]
            axx.hist(xdata, range=(0, xmax), color=color, bins=bins)
            axy.hist(ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins)

            ax.scatter(xdata, ydata, color=color, s=0.2, alpha=0.5)
            # ax.set_xlim(0, xmax)
    
    plot('nde', 'energies')
    if 'dip_trans' in inter_state:
        plot('ntd', 'dip_trans')

    return axs

# %%

def mol_to_png(mol, width=320, height=240):
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, height)
    
    d.drawOptions().setBackgroundColour((1, 1, 1, 0))
    d.drawOptions().padding = 0.05

    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()

def show_atXYZ(atXYZ, name='', smiles=None, inchi=None, skeletal=True, ax=None):
    fig, ax = pca_biplot.figax(ax)
    
    mol = pca_biplot.xyz_to_mol(atXYZ)
    smol = rdkit.Chem.RemoveHs(mol)
    rdkit.Chem.RemoveStereochemistry(smol)
    smiles = rdkit.Chem.MolToSmiles(smol) if smiles is None else smiles
    inchi = rdkit.Chem.MolToInchi(smol) if inchi is None else inchi
    
    png = mol_to_png(rdkit.Chem.RemoveHs(mol) if skeletal else mol)
    pca_biplot.mpl_imshow_png(ax, png)
    ax.set_title(name)
    ax.axis('on')
    ax.get_yaxis().set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_xlabel(f"SMILES={smiles}\n{inchi}", wrap=True)
    print(smiles, inchi)
    # axy.tick_params(axis="y", labelleft=False)
    return ax

# %%

def plot_datasheet(
  per_state, inter_state,
  pops, sgroups,
  noodle, hops,
  delta_E, fosc_time,
  atXYZ, name, smiles=None, inchi=None,
  axs=None, fig=None,
  include_hist=False,
  skeletal=True):
    if axs is None:
        mosaic = [
          ['sg', 'pca', 'pca'],
          ['t0', 'pca', 'pca'],
          ['t1', 'mol', 'pop'],
          ['se', 'nde', 'de'],
          ['t2', 'ntd', 'ft']],
        if include_hist:
            mosaic = [['energies', 'forces', 'dip_perm']] + mosaic
        fig, axs = plt.subplot_mosaic(mosaic, layout='constrained')
    else:
        fig = fig or list(axs.values())[0].figure

    col_state = ['#4DAD15', '#AD2915', '#7515AD']
    dcol_state = {s: c for s, c in zip(pops.state.values, col_state)}
    col_inter = ['#2c3e50', '#C4A000', '#7E5273']
    dcol_inter = {sc: c for sc, c in zip(inter_state.statecomb.values, col_inter)}

    # axs['fde'].sharex(axs['nde'])
    # for i in range(3):
    #     axs[f't{i}'].sharex(axs['nde'])

    if include_hist:
        plot_per_state_histograms(per_state, dcol_state, axs)
    plot_timeplots(pops, delta_E, fosc_time, dcol_state, dcol_inter, axs)
    if sgroups is not None:
        subset = {k: v for k, v in axs.items() if k in ['sg','t0','t1','se','t2','cb_hist','cb_spec']}
        plot_separated_spectra_and_hists(inter_state, sgroups, dcol_inter, axs=subset)
    pca_biplot.plot_noodleplot(noodle, hops, axs['pca'])
    plot_nacs_histograms(inter_state, hops.frame, col_inter, axs)
    show_atXYZ(atXYZ, name, smiles, inchi, skeletal, axs['mol'])
    
    vscale = 1 if include_hist else 5/6
    fig.set_size_inches(8.27,11.69*vscale) # portrait A4
    fig.set_dpi(600)
    plt.show()

    return fig, axs

def create_subfigures(include_hist=False, borders=False):
    def f(sfname, sgs, nrows, ncols, *axnames, **kws):
        nonlocal sfs, axs
        sfs[sfname] = fig.add_subfigure(sgs)
        sfs[sfname].set_facecolor('w')
        axlist = sfs[sfname].subplots(nrows, ncols, **kws)
        for n, ax in zip(axnames, np.atleast_1d(axlist)):
            axs[n] = ax
    
    nrows = 6 if include_hist else 5
    s = 1 if include_hist else 0

    fig, oaxs = plt.subplots(nrows, 3, layout='constrained')
    vscale = 1 if include_hist else 5/6
    fig.set_size_inches(8.27,11.69*vscale) # portrait A4
    if borders:
        fig.set_facecolor('#ddd')
    gs = oaxs[0,0].get_subplotspec().get_gridspec()
    for ax in oaxs.ravel():
        ax.remove()
    
    sfs = {}; axs = {}
    if include_hist:
        f('hist', gs[0,:], 1, 3, 'energies', 'forces', 'dip_perm')
    f('time', gs[s+2:,2], 3, 1, 'pop', 'de', 'ft')
    f('pca', gs[s+0:s+2,1:], 1, 1, 'pca')
    f('de', gs[s+0:,0], 5+2, 1,
      'sg', 't0', 't1', 'se', 't2', 'cb_spec', 'cb_hist',
      height_ratios=([1]*5)+([0.1]*2))
    f('nacs', gs[s+3:,1], 2, 1, 'nde', 'ntd')
    f('mol', gs[s+2,1], 1, 1, 'mol')
    return fig, sfs, axs