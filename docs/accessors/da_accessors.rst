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

    da.st.to_xyz()

The :py:class:`shnitsel.xarray.DAShnitselAccessor` is opinionated about whether a method
is suitable for the :py:class:`xarray.DataArray` in question.  The list of suitable
When working in a notebook, the HTML representation of the accessor (type ``da.sh``) will show
suitable methods and explain why the others are unsuitable. This information is also
available through the ``da.st.suitable`` list and the ``da.st.unsuitable`` dictionary.

.. currentmodule:: xarray
.. autosummary::
    :toctree: ../_generated
    :template: autosummary/accessor_method.rst

    DataArray.st.angle
    DataArray.st.assign_levels
    DataArray.st.calc_ci
    DataArray.st.default_mol
    DataArray.st.dihedral
    DataArray.st.distance
    DataArray.st.expand_midx
    DataArray.st.flatten_levels
    DataArray.st.get_bond_lengths
    DataArray.st.get_hop_types
    DataArray.st.hop_indices
    DataArray.st.keep_norming
    DataArray.st.last_time_where
    DataArray.st.mgroupby
    DataArray.st.msel
    DataArray.st.norm
    DataArray.st.pairwise_dists_pca
    DataArray.st.pca
    DataArray.st.relativize
    DataArray.st.sel_trajs
    DataArray.st.subtract_combinations
    DataArray.st.sudi
    DataArray.st.time_grouped_ci
    DataArray.st.to_mol
    DataArray.st.to_xyz
    DataArray.st.traj_to_xyz
    DataArray.st.trajs_with_hops
    DataArray.st.ts_to_time