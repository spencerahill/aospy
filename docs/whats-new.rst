What's New
==========

v0.1.1
------
This release includes fixes for a number of bugs mistakenly introduced in the
refactoring of the variable loading step of ``calc.py`` (:issue:`90`), as well as
support for xarray version 0.9.0.

Ehancements
~~~~~~~~~~~
- Support for xarray version 0.9.0.  By `Spencer Clark <https://github.com/spencerkclark>`_.

Bug fixes
~~~~~~~~~
- Fix an instance where the name for pressure half levels was mistakenly
  replaced with the name for the pressure full levels (:issue:`126`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Prevent workaround for dates outside the ``pd.Timestamp`` valid range from
  being applied to dates within the ``pd.Timestamp`` valid range (:issue:`128`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Ensure that all DataArrays associated with :py:class:`aospy.Var` objects have a time
  weights coordinate with CF-compliant time units.  This allows them to be cast
  as the type ``np.timedelta64``, and be safely converted to have units of days before
  taking time-weighted averages (:issue:`128`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix a bug where the time weights were not subset in time prior to taking a time weighted average; this caused computed
  seasonal averages to be too small.  To prevent this from failing silently again, 
  we now raise a ``ValueError`` if the time coordinate of the time weights
  is not identical to the time coordinate of the array associated with the
  :py:class:`aospy.Var` (:issue:`128`).  By `Spencer Clark <https://github.com/spencerkclark>`_.

v0.1 (24 January 2017)
----------------------
- Initial release!
- Contributors:
  
  - `Spencer Hill <https://github.com/spencerahill>`_
  - `Spencer Clark <https://github.com/spencerkclark>`_
