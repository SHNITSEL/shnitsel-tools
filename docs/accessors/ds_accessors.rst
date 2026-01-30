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

    ds.st.to_xyz()

The :py:class:`shnitsel.xarray.DSShnitselAccessor` is opinionated about whether a method
is suitable for the :py:class:`xarray.Dataset` in question.  The list of suitable
When working in a notebook, the HTML representation of the accessor (type ``ds.sh``) will show
suitable methods and explain why the others are unsuitable. This information is also
available through the ``ds.st.suitable`` list and the ``ds.st.unsuitable`` dictionary.

.. currentmodule:: xarray
.. autosummary::
    :toctree: ../_generated
    :template: autosummary/accessor_method.rst

    Dataset.st.assign_fosc
    Dataset.st.assign_levels
    Dataset.st.ds_broaden_gauss
    Dataset.st.calc_classical_populations
    Dataset.st.default_mol
    Dataset.st.calculate_energy_filtranda
    Dataset.st.expand_midx
    Dataset.st.find_hops
    Dataset.st.flatten_levels
    Dataset.st.get_cutoffs
    Dataset.st.get_inter_state
    Dataset.st.get_per_state
    Dataset.st.iconds_to_frames
    Dataset.st.mgroupby
    Dataset.st.msel
    Dataset.st.pca_and_hops
    Dataset.st.save_frames
    Dataset.st.sel_trajs
    Dataset.st.setup_frames
    Dataset.st.spectra_all_times
    Dataset.st.truncate
    Dataset.st.ts_to_time
    Dataset.st.validate
    Dataset.st.write_ase_db

..
    missing
    Dataset.st.assign_fosc
    Dataset.st.ds_broaden_gauss
    Dataset.st.find_hops
    Dataset.st.get_cutoffs
    Dataset.st.iconds_to_frames
    Dataset.st.save_frames
    Dataset.st.setup_frames
    Dataset.st.spectra_all_times
    Dataset.st.ts_to_time