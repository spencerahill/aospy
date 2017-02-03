What's New
==========

v0.1.1
------
This release includes fixes for a number of bugs mistakenly introduced in the
refactoring of the variable loading step of ``calc.py`` (:pull:`90`), as well as
support for xarray version 0.9.1.

Ehancements
~~~~~~~~~~~
- Support for xarray version 0.9.1 and require it or a later xarray version.  By `Spencer Clark <https://github.com/spencerkclark>`_ and `Spencer Hill <https://github.com/spencerahill>`_.
- Better support for variable names relating to "bounds" dimension of input data files.  "bnds", "bounds", and "nv" now all supported (:pull:`140`).  By `Spencer Hill <https://github.com/spencerahill>`_.
- When coercing dims of input data to aospy's internal names, for scalars change only the name; for non-scalars change the name, force them to have a coord, and copy over their attrs (:pull:`140`).  By `Spencer Hill <https://github.com/spencerahill>`_.

Bug fixes
~~~~~~~~~
- Fix bug involving loading data that has dims that lack coords (which is possible as of xarray v0.9.0).  By `Spencer Hill <https://github.com/spencerahill>`_.
- Fix an instance where the name for pressure half levels was mistakenly
  replaced with the name for the pressure full levels (:pull:`126`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Prevent workaround for dates outside the ``pd.Timestamp`` valid range from
  being applied to dates within the ``pd.Timestamp`` valid range (:pull:`128`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Ensure that all DataArrays associated with :py:class:`aospy.Var` objects have a time
  weights coordinate with CF-compliant time units.  This allows them to be cast
  as the type ``np.timedelta64``, and be safely converted to have units of days before
  taking time-weighted averages (:pull:`128`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix a bug where the time weights were not subset in time prior to taking a time weighted average; this caused computed
  seasonal averages to be too small.  To prevent this from failing silently again, 
  we now raise a ``ValueError`` if the time coordinate of the time weights
  is not identical to the time coordinate of the array associated with the
  :py:class:`aospy.Var` (:pull:`128`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Enable calculations to be completed using data saved as a single time-slice
  on disk (fixes :issue:`132` through :pull:`135`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix bug where workaround for dates outside the ``pd.Timestamp`` valid range
  caused a mismatch between the data loaded and the data requested (fixes
  :issue:`138` through :pull:`139`). By `Spencer Clark <https://github.com/spencerkclark>`_.

v0.1 (24 January 2017)
----------------------
- Initial release!
- Contributors:
  
  - `Spencer Hill <https://github.com/spencerahill>`_
  - `Spencer Clark <https://github.com/spencerkclark>`_
