.. currentmodule:: shnitsel

DataArray accessor methods
==========================

Each of the following functions can be used in an object-oriented style thanks
to :py:mod:`xarray`'s accessor mechanism. The class :py:class:`shnitsel.xarray.DAShnitselAccessor`
is registered as an accessor upon import of :py:mod:`shnitsel.xarray`, which means
that it is instatiated as the ``sh`` property of every :py:class:`xarray.DataArray`.

In practice, this means that instead of writing

.. code-block::

    shnitsel.core.postprocess.to_xyz(da)

you can write

.. code-block::

    da.sh.to_xyz()

The :py:class:`shnitsel.xarray.DAShnitselAccessor` is opinionated about whether functions
are relevant for the :py:class:`xarray.DataArray` in question. Irrelevant functions
are made unavailable, in the sense that they do not appear in ``__dir__()`` and 
calling them will raise a :py:class:`TypeError`.
In most cases, unavailable functions would have failed anyway. When working in a
notebook, the HTML representation of the accessor (type ``da.sh``) will show
available methods and explain why the others are unavailable.

.. autosummary::
    :toctree: ../_generated

    shnitsel.core.postprocess.angle
    shnitsel.core.xrhelpers.assign_levels
    shnitsel.core.postprocess.calc_ci
    shnitsel.core.postprocess.default_mol
    shnitsel.core.postprocess.dihedral
    shnitsel.core.postprocess.distance
    shnitsel.core.xrhelpers.expand_midx
    shnitsel.core.xrhelpers.flatten_levels
    shnitsel.core.filtre.get_bond_lengths
    shnitsel.core.postprocess.get_hop_types
    shnitsel.core.postprocess.hop_indices
    shnitsel.core.postprocess.keep_norming
    shnitsel.core.filtre.last_time_where
    shnitsel.core.xrhelpers.mgroupby
    shnitsel.core.xrhelpers.msel
    shnitsel.core.postprocess.norm
    shnitsel.core.postprocess.pairwise_dists_pca
    shnitsel.core.postprocess.pca
    shnitsel.core.postprocess.relativize
    shnitsel.core.xrhelpers.sel_trajs
    shnitsel.core.postprocess.subtract_combinations
    shnitsel.core.postprocess.sudi
    shnitsel.core.postprocess.time_grouped_ci
    shnitsel.core.postprocess.to_mol
    shnitsel.core.postprocess.to_xyz
    shnitsel.core.postprocess.traj_to_xyz
    shnitsel.core.postprocess.trajs_with_hops
    shnitsel.core.postprocess.ts_to_time