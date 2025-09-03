.. currentmodule:: shnitsel

Dataset accessor methods
========================


Each of the following functions can be used in an object-oriented style thanks
to :py:mod:`xarray`'s accessor mechanism. The class :py:class:`shnitsel.xarray.DSShnitselAccessor`
is registered as an accessor upon import of :py:mod:`shnitsel.xarray`, which means
that it is instatiated as the ``sh`` property of every :py:class:`xarray.Dataset`.

In practice, this means that instead of writing

.. code-block::

    shnitsel.core.postprocess.to_xyz(ds)

you can write

.. code-block::

    ds.sh.to_xyz()

The :py:class:`shnitsel.xarray.DSShnitselAccessor` is opinionated about whether functions
are relevant for the :py:class:`xarray.Dataset` in question. Irrelevant functions
are made unavailable, in the sense that they do not appear in ``__dir__()`` and 
calling them will raise a :py:class:`TypeError`.
In most cases, unavailable functions would have failed anyway. When working in a
notebook, the HTML representation of the accessor (type ``ds.sh``) will show
available methods and explain why the others are unavailable.

.. autosummary::
    :toctree: ../_generated

    shnitsel.core.postprocess.assign_fosc
    shnitsel.core.xrhelpers.assign_levels
    shnitsel.core.postprocess.broaden_gauss
    shnitsel.core.postprocess.calc_pops
    shnitsel.core.postprocess.default_mol
    shnitsel.core.filtre.energy_filtranda
    shnitsel.core.xrhelpers.expand_midx
    shnitsel.core.postprocess.find_hops
    shnitsel.core.xrhelpers.flatten_levels
    shnitsel.core.filtre.get_cutoffs
    shnitsel.core.postprocess.get_inter_state
    shnitsel.core.postprocess.get_per_state
    shnitsel.core.parse.sharc_icond.iconds_to_frames
    shnitsel.core.xrhelpers.mgroupby
    shnitsel.core.xrhelpers.msel
    shnitsel.core.postprocess.pca_and_hops
    shnitsel.core.xrhelpers.save_frames
    shnitsel.core.xrhelpers.sel_trajs
    shnitsel.core.postprocess.setup_frames
    shnitsel.core.plot.spectra3d.spectra_all_times
    shnitsel.core.filtre.truncate
    shnitsel.core.postprocess.ts_to_time
    shnitsel.core.postprocess.validate
    shnitsel.core.ase.write_ase
