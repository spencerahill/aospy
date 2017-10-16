.. _whats-new:

What's New
==========

.. _whats-new.0.2.1:

v0.2.1 (unreleased)
-------------------

Documentation
~~~~~~~~~~~~~

Corrected link to documentation badge on repository main page (:pull:213). By DaCoEx <https://github.com/dacoex>_.

=======
Enhancements
~~~~~~~~~~~~

- Remove potentially confusing attributes from example netcdf files.
  (closes :issue:`214` via :pull:`216`). By `Micah Kim
  <https://github.com/micahkim23>`_.

Dependencies
~~~~~~~~~~~~

- ``aospy`` now requires a minimum version of ``distributed`` of
  1.17.1 (fixes :issue:`210` via :pull:`211`).

.. _whats-new.0.2:

v0.2 (26 September 2017)
------------------------

This release includes some new features plus several bugfixes.  The
bugfixes include some that previously made using aospy on
pressure-interpolated data very problematic.  We have also improved
support for reading in data from the WRF and CAM atmospheric models.

As of this release, aospy has at least 2(!) confirmed regular users
that aren't the original aospy developers, bringing the worldwide
total of users up to at least 4.  The first user-generated Github
Issues have now also been created.  We're a real thing!

Enhancements
~~~~~~~~~~~~

- Use ``dask.bag`` coupled with ``dask.distributed`` rather than
  ``multiprocess`` to parallelize computations (closes :issue:`169`
  via :pull:`172`).  This enables the optional use of an external
  ``distributed.Client`` to leverage computational resources across
  multiple nodes of a cluster. By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Improve support for WRF and NCAR CAM model data by adding the
  internal names they use for grid attributes to aospy's lists of
  potential names to search for.  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- Allow a user to specify a custom preprocessing function in all
  DataLoaders to prepare data for processing with aospy.  This could
  be used, for example, to add a CF-compliant units attribute to the
  time coordinate if it is not present in a set of files.  Addresses
  :issue:`177` via :pull:`180`.  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Remove ``dask.async`` import in ``model.py``; no longer needed, and
  also prevents warning message from dask regarding location of
  ``get_sync`` function  (:pull:`195`).  By
  `Spencer Hill <https://github.com/spencerahill>`_.


Dependencies
~~~~~~~~~~~~

- ``multiprocess`` is no longer required for submitting ``aospy`` calculations
  in parallel (see discussion in :issue:`169` and pull request :pull:`172`).
- ``aospy`` now requires an installation of ``dask`` with version
  greater than or equal to 0.14 (see discussion in pull request
  :pull:`172`).

Bug Fixes
~~~~~~~~~

- Remove faulty logic for calculations with data coming from multiple
  runs.  Eventually this feature will be properly implemented (fixes
  :issue:`117` via :pull:`178`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- Only run tests that require optional dependencies if those
  dependencies are actually installed (fixes :issue:`167` via
  :pull:`176`).  By `Spencer Hill <https://github.com/spencerahill>`_.
- Remove obsolete ``operator.py`` module (fixes :issue:`174` via
  :pull:`175`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Fix workaround for dates with years less than 1678 to support units
  attributes with a reference date years not equal to 0001 (fixes
  :issue:`188` via :pull:`189`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix bug which would prevent users from analyzing a subset within the
  Timestamp-valid range from a dataset which
  included data from outside the Timestamp-valid range (fixed in
  :pull:`189`). By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Toggle the ``mask_and_scale`` option to ``True`` when reading in netCDF files
  to enable missing values encoded as floats to be converted to NaN's
  (fixes :issue:`190` via :pull:`192`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Force regional calculations to mask gridcell weights where the loaded
  datapoints were invalid instead of just masking points outside the desired
  region (fixes :issue:`190` via :pull:`192`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Retain original input data's mask during gridpoint-by-gridpoint
  temporal averages (fixes :issue:`193` via :pull:`196`).  By `Spencer
  Hill <https://github.com/spencerahill>`_.
- Always write output to a tar file in serial to prevent empty header file
  errors (fixes :issue:`75` via :pull:`197`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Allow ``aospy`` to use grid attributes that are only defined in ``Run``
  objects. Previously if a grid attribute were defined only in a ``Run``
  object and not also in the Run's corresponding ``Model``, an error would
  be raised (fixes :issue:`187` via :pull:`199`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- When input data for a calculation has a time bounds array, overwrite
  its time array with the average of the start and end times for each
  timestep.  Prevents bug wherein time arrays equal to either the
  start or end bounds get mistakenly grouped into the wrong time
  interval, i.e. the wrong month or year (fixes :issue `185` via
  :pull:`200`).  By `Spencer Hill <https://github.com/spencerahill>`_.

.. _whats-new.0.1.2:

v0.1.2 (30 March 2017)
----------------------

This release improves the process of submitting multiple calculations
for automatic execution.  The user interface, documentation, internal
logic, and packaging all received upgrades and/or bugfixes.

We also now have a `mailing list`_.  Join it to follow and/or post
your own usage questions, bug reports, suggestions, etc.

.. _mailing list: https://groups.google.com/d/forum/aospy

Enhancements
~~~~~~~~~~~~

- Include an example library of aospy objects that works
  out-of-the-box with the provided example main script (:pull:`155`).
  By `Spencer Clark <https://github.com/spencerkclark>`_ and `Spencer
  Hill <https://github.com/spencerahill>`_.
- Improve :ref:`examples` page of the documentation by using this new
  example object library (:pull:`164`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- Improve readability/usability of the included example script
  ``aospy_main.py`` for submitting aospy calculations by moving all
  internal logic into new ``automate.py`` module (:pull:`155`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_ and `Spencer
  Hill <https://github.com/spencerahill>`_.
- Enable user to specify whether or not to write output to .tar files
  (in addition to the standard output).  Also document an error that
  occurs when writing output to .tar files for sufficiently old
  versions of tar (including the version that ships standard on
  MacOS), and print a warning when errors are caught during the 'tar'
  call (:pull:`160`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.

Bug fixes
~~~~~~~~~

- Update packaging specifications such that the example main script
  and tutorial notebook actually ship with aospy as intended (fixes
  :issue:`149` via :pull:`161`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- Use the 'scipy' engine for the `xarray.DataArray.to_netcdf`_
  call when writing aospy calculation outputs to disk to prevent a bug
  when trying to re-write to an existing netCDF file (fixes
  :issue:`157` via :pull:`160`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.

.. _xarray.DataArray.to_netcdf : http://xarray.pydata.org/en/stable/generated/xarray.DataArray.to_netcdf.html

.. _whats-new.0.1.1:

v0.1.1 (2 March 2017)
---------------------

This release includes fixes for a number of bugs mistakenly introduced
in the refactoring of the variable loading step of ``calc.py``
(:pull:`90`), as well as support for xarray version 0.9.1.

Enhancements
~~~~~~~~~~~~
- Support for xarray version 0.9.1 and require it or a later xarray
  version.  By `Spencer Clark <https://github.com/spencerkclark>`_ and
  `Spencer Hill <https://github.com/spencerahill>`_.
- Better support for variable names relating to "bounds" dimension of
  input data files.  "bnds", "bounds", and "nv" now all supported
  (:pull:`140`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- When coercing dims of input data to aospy's internal names, for
  scalars change only the name; for non-scalars change the name, force
  them to have a coord, and copy over their attrs (:pull:`140`).  By
  `Spencer Hill <https://github.com/spencerahill>`_.

Bug fixes
~~~~~~~~~
- Fix bug involving loading data that has dims that lack coords (which
  is possible as of xarray v0.9.0).  By `Spencer Hill
  <https://github.com/spencerahill>`_.
- Fix an instance where the name for pressure half levels was
  mistakenly replaced with the name for the pressure full levels
  (:pull:`126`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Prevent workaround for dates outside the ``pd.Timestamp`` valid
  range from being applied to dates within the ``pd.Timestamp`` valid
  range (:pull:`128`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Ensure that all DataArrays associated with :py:class:`aospy.Var`
  objects have a time weights coordinate with CF-compliant time units.
  This allows them to be cast as the type ``np.timedelta64``, and be
  safely converted to have units of days before taking time-weighted
  averages (:pull:`128`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Fix a bug where the time weights were not subset in time prior to
  taking a time weighted average; this caused computed seasonal
  averages to be too small.  To prevent this from failing silently
  again, we now raise a ``ValueError`` if the time coordinate of the
  time weights is not identical to the time coordinate of the array
  associated with the :py:class:`aospy.Var` (:pull:`128`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Enable calculations to be completed using data saved as a single
  time-slice on disk (fixes :issue:`132` through :pull:`135`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix bug where workaround for dates outside the ``pd.Timestamp``
  valid range caused a mismatch between the data loaded and the data
  requested (fixes :issue:`138` through :pull:`139`). By `Spencer
  Clark <https://github.com/spencerkclark>`_.

.. _whats-new.0.1:

v0.1 (24 January 2017)
----------------------
- Initial release!
- Contributors:

  - `Spencer Hill <https://github.com/spencerahill>`_
  - `Spencer Clark <https://github.com/spencerkclark>`_
