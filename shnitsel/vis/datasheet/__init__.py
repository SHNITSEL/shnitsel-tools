# noqa: F401

# from typing import TYPE_CHECKING
# import matplotlib.pyplot as plt
# import numpy as np

# if TYPE_CHECKING:
#     from matplotlib.figure import SubFigure
#     from matplotlib.axes import Axes

# from ..plot import pca_biplot as pca_biplot

# from ..spectra import (
#     get_spectrum as get_spectrum,
#     calc_spectra as calc_spectra,
#     get_spectra_groups as get_spectra_groups,
#     sep_ground_excited_spectra as sep_ground_excited_spectra,
# )

# from .per_state_hist import plot_per_state_histograms as plot_per_state_histograms
# from .time import plot_pops as plot_pops, plot_timeplots as plot_timeplots
# from .dip_trans_hist import (
#     single_hist as single_hist,
#     plot_dip_trans_histograms as plot_dip_trans_histograms,
#     plot_spectra as plot_spectra,
#     plot_separated_spectra_and_hists as plot_separated_spectra_and_hists,
#     plot_separated_spectra_and_hists_groundstate as plot_separated_spectra_and_hists_groundstate,
# )
# from .nacs_hist import plot_nacs_histograms as plot_nacs_histograms


from .datasheet import Datasheet as Datasheet

__all__ = ['Datasheet']  # , 'show_atXYZ']


# def plot_datasheet(
#     per_state,
#     inter_state,
#     pops,
#     sgroups,
#     noodle,
#     hops,
#     delta_E,
#     fosc_time,
#     atXYZ,
#     name,
#     smiles=None,
#     inchi=None,
#     axs=None,
#     fig=None,
#     include_hist=False,
#     skeletal=True,
# ):
#     if axs is None:
#         mosaic = [
#             ('sg', 'pca', 'pca'),
#             ('t0', 'pca', 'pca'),
#             ('t1', 'mol', 'pop'),
#             ('se', 'nde', 'de'),
#             ('t2', 'ntd', 'ft'),
#         ]
#         if include_hist:
#             mosaic = [('energy', 'forces', 'dip_perm')] + mosaic
#         fig, axs = plt.subplot_mosaic(mosaic, layout='constrained')  # type: ignore
#     else:
#         fig = fig or list(axs.values())[0].figure

#     col_state = ['#4DAD15', '#AD2915', '#7515AD']
#     dcol_state = {s: c for s, c in zip(pops.state.values, col_state)}
#     col_inter = ['#2c3e50', '#C4A000', '#7E5273']
#     dcol_inter = {sc: c for sc, c in zip(inter_state.statecomb.values, col_inter)}

#     # axs['fde'].sharex(axs['nde'])
#     # for i in range(3):
#     #     axs[f't{i}'].sharex(axs['nde'])

#     if include_hist:
#         plot_per_state_histograms(per_state, dcol_state, axs)
#     plot_timeplots(pops, delta_E, fosc_time, dcol_state, dcol_inter, axs)
#     if sgroups is not None:
#         subset = {
#             k: v
#             for k, v in axs.items()
#             if k in ['sg', 't0', 't1', 'se', 't2', 'cb_hist', 'cb_spec']
#         }
#         plot_separated_spectra_and_hists(inter_state, sgroups, dcol_inter, axs=subset)
#     pca_biplot.plot_noodleplot(noodle, hops, axs['pca'])
#     plot_nacs_histograms(inter_state, hops.frame, col_inter, axs)
#     show_atXYZ(atXYZ, name, smiles, inchi, skeletal, axs['mol'])

#     vscale = 1 if include_hist else 5 / 6
#     fig.set_size_inches(8.27, 11.69 * vscale)  # portrait A4
#     fig.set_dpi(600)
#     plt.show()

#     return fig, axs


# def create_subfigures(include_hist=False, borders=False):
#     sfs: dict[str, SubFigure] = {}
#     axs: dict[str, Axes] = {}

#     def f(sfname, sgs, nrows, ncols, *axnames, **kws):
#         nonlocal sfs, axs
#         sfs[sfname] = fig.add_subfigure(sgs)
#         sfs[sfname].set_facecolor('w')
#         axlist = sfs[sfname].subplots(nrows, ncols, **kws)
#         for n, ax in zip(axnames, np.atleast_1d(axlist)):
#             axs[n] = ax

#     nrows = 6 if include_hist else 5
#     s = 1 if include_hist else 0

#     fig, oaxs = plt.subplots(nrows, 3, layout='constrained')
#     vscale = 1 if include_hist else 5 / 6
#     fig.set_size_inches(8.27, 11.69 * vscale)  # portrait A4
#     if borders:
#         fig.set_facecolor('#ddd')
#     gs = oaxs[0, 0].get_subplotspec().get_gridspec()
#     for ax in oaxs.ravel():
#         ax.remove()

#     if include_hist:
#         f('hist', gs[0, :], 1, 3, 'energy', 'forces', 'dip_perm')
#     f('time', gs[s + 2 :, 2], 3, 1, 'pop', 'de', 'ft')
#     f('pca', gs[s + 0 : s + 2, 1:], 1, 1, 'pca')
#     hr = ([1] * 5) + ([0.1] * 2)
#     f('de', gs[s + 0 :, 0], 5 + 2, 1,
#       'sg', 't0', 't1', 'se', 't2', 'cb_spec', 'cb_hist',
#        height_ratios=hr)  # fmt: skip
#     f('nacs', gs[s + 3 :, 1], 2, 1, 'nde', 'ntd')
#     f('mol', gs[s + 2, 1], 1, 1, 'mol')
#     return fig, sfs, axs
