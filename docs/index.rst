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

Reading data from...
--------------------

...trajectories:

.. code-block::

      # One of the following:
      ds = st.read_trajs('path/to/Singlet_1/', kind='sharc')
      ds = st.read_trajs('path/to/Newton-X/', kind='nx')
      ds = st.read_trajs('path/to/PyRAI2MD/', kind='pyrai2md')

      # And then:
      ds = ds.sh.setup_frames()
   
...ASE databases:

.. code-block::

      # One of the following:
      ds = st.read_ase('path/to/ase.db', kind='spainn')
      ds = st.read_ase('path/to/ase.db', kind='schnet')
      ds = st.read_ase('path/to/geometries_only.db', kind=None)

...SHNITSEL-style NetCDF4:

.. code-block::

      ds = st.open_frames('path/to/dataset.nc')


Saving data to...
-----------------

...NetCDF4:

.. code-block::

   ds.sh.save_frames('path/to/dataset.nc')

...ASE:

.. code-block::

      ds = st.write_ase('path/to/ase.db', kind='spainn')
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
   priority
   api/index