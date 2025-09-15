Dataset accessor methods
========================


Each of the following functions can be used in an object-oriented style thanks
to :py:mod:`xarray`'s accessor mechanism. The class :py:class:`shnitsel._generated_accessors.DatasetAccessor`
is registered as an accessor upon import of :py:mod:`shnitsel.xarray`, which means
that it is instatiated as the ``sh`` property of every :py:class:`xarray.Dataset`.

In practice, this means that instead of writing

.. code-block::

    shnitsel.core.postprocess.to_xyz(ds)

you can write

.. code-block::

    ds.sh.to_xyz()

The :py:class:`shnitsel.xarray.DSShnitselAccessor` is opinionated about whether a method
is suitable for the :py:class:`xarray.Dataset` in question.  The list of suitable
When working in a notebook, the HTML representation of the accessor (type ``ds.sh``) will show
suitable methods and explain why the others are unsuitable. This information is also
available through the ``ds.sh.suitable`` list and the ``ds.sh.unsuitable`` dictionary.

.. currentmodule:: xarray
.. autosummary::
    :toctree: ../_generated
    :template: autosummary/accessor_method.rst

    Dataset.sh.assign_fosc
    Dataset.sh.assign_levels
    Dataset.sh.ds_broaden_gauss
    Dataset.sh.calc_pops
    Dataset.sh.default_mol
    Dataset.sh.energy_filtranda
    Dataset.sh.expand_midx
    Dataset.sh.find_hops
    Dataset.sh.flatten_levels
    Dataset.sh.get_cutoffs
    Dataset.sh.get_inter_state
    Dataset.sh.get_per_state
    Dataset.sh.iconds_to_frames
    Dataset.sh.mgroupby
    Dataset.sh.msel
    Dataset.sh.pca_and_hops
    Dataset.sh.save_frames
    Dataset.sh.sel_trajs
    Dataset.sh.setup_frames
    Dataset.sh.spectra_all_times
    Dataset.sh.truncate
    Dataset.sh.ts_to_time
    Dataset.sh.validate
    Dataset.sh.write_ase
