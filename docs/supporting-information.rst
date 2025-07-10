Functions mentioned in supporting information
=============================================

The following are mentioned in the SI:

.. autofunction:: shnitsel.parse.sharc_icond.dirs_of_iconds

----

.. autofunction:: shnitsel.parse.read_trajs

The same function is exposed as ``shnitsel.read_trajs``,
and it might be better to change the SI to use this,
at which point we hardly need to expose ``shnitsel.parse``.
Oh, now I think about it, there should have been just a single
function called ``read`` or ``parse`` that can handle anything
based on the ``kind`` parameter...

----

.. autofunction:: shnitsel.open_frames

.. .. automethod:: shnitsel.xarray.DAShnitselAccessor.get_bond_lengths

----

We'll have to refactor the accessor methods to document them properly.

But for now:

.. automethod:: shnitsel.core.filtre.get_bond_lengths


.. automethod:: shnitsel.core.xrhelpers.msel


.. automethod:: shnitsel.core.xrhelpers.sel_trajs


.. automethod:: shnitsel.core.xrhelpers.save_frames


.. automethod:: shnitsel.core.postprocess.dihedral

.. automethod:: shnitsel.core.postprocess.get_inter_state

.. automethod:: shnitsel.core.plot.spectra3d.spectra_all_times

----

And now some bespoke plotting functions.

.. autofunction:: shnitsel.plot.biplot_kde

.. autofunction:: shnitsel.plot.ski_plots

.. autoclass:: shnitsel.plot.Datasheet
    :members: 

----

Next the SI mentions  

.. .. autofunction:: shnitsel.