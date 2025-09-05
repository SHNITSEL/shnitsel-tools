.. currentmodule:: shnitsel

DataArray accessor methods
==========================

Each of the following functions can be used in an object-oriented style thanks
to :py:mod:`xarray`'s accessor mechanism. The class :py:class:`shnitsel._generated_accessors.DataArrayAccessor`
is registered as an accessor upon import of :py:mod:`shnitsel.xarray`, which means
that it is instatiated as the ``sh`` property of every :py:class:`xarray.DataArray`.

In practice, this means that instead of writing

.. code-block::

    shnitsel.core.postprocess.to_xyz(da)

you can write

.. code-block::

    da.sh.to_xyz()

The :py:class:`shnitsel.xarray.DAShnitselAccessor` is opinionated about whether a method
is suitable for the :py:class:`xarray.DataArray` in question.  The list of suitable
When working in a notebook, the HTML representation of the accessor (type ``da.sh``) will show
suitable methods and explain why the others are unsuitable. This information is also
available through the ``da.sh.suitable`` list and the ``da.sh.unsuitable`` dictionary.

.. currentmodule:: xarray
.. autosummary::
    :toctree: ../_generated
    :template: autosummary/accessor_method.rst

    DataArray.sh.angle
    DataArray.sh.assign_levels
    DataArray.sh.calc_ci
    DataArray.sh.default_mol
    DataArray.sh.dihedral
    DataArray.sh.distance
    DataArray.sh.expand_midx
    DataArray.sh.flatten_levels
    DataArray.sh.get_bond_lengths
    DataArray.sh.get_hop_types
    DataArray.sh.hop_indices
    DataArray.sh.keep_norming
    DataArray.sh.last_time_where
    DataArray.sh.mgroupby
    DataArray.sh.msel
    DataArray.sh.norm
    DataArray.sh.pairwise_dists_pca
    DataArray.sh.pca
    DataArray.sh.relativize
    DataArray.sh.sel_trajs
    DataArray.sh.subtract_combinations
    DataArray.sh.sudi
    DataArray.sh.time_grouped_ci
    DataArray.sh.to_mol
    DataArray.sh.to_xyz
    DataArray.sh.traj_to_xyz
    DataArray.sh.trajs_with_hops
    DataArray.sh.ts_to_time