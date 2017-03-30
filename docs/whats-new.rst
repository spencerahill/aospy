.. _whats-new:

What's New
==========

.. _whats-new.0.1.2:

v0.1.2 (XX March 2017)
----------------------

This release improves the process of submitting multiple calculations
for automatic execution.  The user interface, documentation, internal
logic, and packaging all received upgrades and/or bugfixes.

We also now have a `mailing list`_.  Join it to follow and/or post
your own usage questions, bug reports, suggestions, etc.

.. _mailing list: https://groups.google.com/d/forum/aospy

Enhancements
~~~~~~~~~~~~

- Improve Examples page of the documentation by using the newly
  improved example main script and object libraries (:pull:`164`).  By
  `Spencer Hill <https://github.com/spencerahill>`_.
- Include an example library of aospy objects that works
  out-of-the-box with the provided example main script (:pull:`155`).
  By `Spencer Clark <https://github.com/spencerkclark>`_ and `Spencer
  Hill <https://github.com/spencerahill>`_.
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
- Use the 'scipy' engine for the `xarray.DataArray.to_netcdf
  <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.to_netcdf.html?highlight=to_netcdf>`_
  call when writing aospy calculation outputs to disk to prevent a bug
  when trying to re-write to an existing netCDF file (fixes
  :issue:`157` via :pull:`160`).  By `Spencer Hill
  <https://github.com/spencerahill>`_.

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
