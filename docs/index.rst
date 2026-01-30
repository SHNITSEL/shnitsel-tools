.. shnitsel-tools documentation master file, created by
   sphinx-quickstart on Mon Jul  7 10:53:00 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

********************************
``shnitsel-tools`` documentation
********************************

Getting started
===============

First install,

.. code-block::

     pip install shnitsel-tools

then import.

.. code-block::

     import shnitsel as st
     import shnitsel.xarray # Activate accessors
     from shnitsel.data.tree import tree_to_frames

Reading data from...
--------------------

...trajectories:

.. code-block::

      # One of the following:
      db = st.read('path/to/Singlet_1/', kind='sharc')
      db = st.read('path/to/Newton-X/', kind='nx')
      db = st.read('path/to/PyRAI2MD/', kind='pyrai2md')
      db = st.read('path/to/ase.db')

      # And then:
      ds = tree_to_frames(db)

...SHNITSEL-style NetCDF4:

.. code-block::

      ds = st.read('path/to/dataset.nc')

.. #COMMENT OUT -- check the following
   ...ASE databases:

   .. code-block::

         # One of the following:
         ds = st.read('path/to/ase.db', kind='spainn')
         ds = st.read('path/to/ase.db', kind='schnet')
         ds = st.read('path/to/geometries_only.db', kind=None)



Saving data to...
-----------------

...NetCDF4:

.. code-block::

   ds.st.write_shnitsel_file('path/to/dataset.nc')


.. #COMMENT OUT -- what is the correct new syntax?
   ...ASE:

   .. code-block::

         ds.sh.write_ase('path/to/ase.db', kind='spainn')
         # see above for other kinds

Select...
---------

...a variable:

.. code-block::

   ds['energy'] # ...or:
   ds.energy

...a trajectory:

.. code-block::

   ds.sel(trajid=1)

...a point in time:

.. code-block::

   ds.sel(time=1)

.. #COMMENT OUT -- TODO
   Derive...
   ---------

   ...bond lengths:
   TODO


Reference
=========

.. toctree::
   :maxdepth: 1

   top-level
   accessors/ds_accessors
   accessors/da_accessors
   api/index
   tutorials