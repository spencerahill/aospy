#############
API Reference
#############

Here we provide the reference documentation for aospy's public API.
If you are new to the package and/or just trying to get a feel for the
overall workflow, you are better off starting in the main
documentation sections.

.. warning::

   aospy is under active development.  While we strive to maintain
   backwards compatibility, it is likely that some breaking changes to
   the codebase will occur in the future as aospy is improved.

Core Hierarchy for Input Data
=============================

aospy provides three classes for specifying the location and
characteristics of data saved on disk as netCDF files that the user
wishes to use as input data for aospy calculations: ``Proj``,
``Model``, and ``Run``.

Proj
----

.. autoclass:: aospy.proj.Proj
    :members:
    :undoc-members:

    .. automethod:: aospy.proj.Proj.__init__

Model
-----

.. autoclass:: aospy.model.Model
    :members:
    :undoc-members:

    .. automethod:: aospy.model.Model.__init__

Run
---

.. autoclass:: aospy.run.Run
    :members:
    :undoc-members:

    .. automethod:: aospy.run.Run.__init__

DataLoaders
===========

``Run`` objects rely on a helper "data loader" to specify how to find
their underlying data that is saved on disk.  This mapping of
variables, time ranges, and potentially other parameters to the
location of the corresponding data on disk can differ among modeling
centers or even between different models at the same center.

Currently supported data loader types are ``DictDataLoader``,
``NestedDictDataLoader``, and ``GFDLDataLoader`` Each of these inherit
from the abstract base ``DataLoader`` class.

.. autoclass:: aospy.data_loader.DataLoader
    :members:
    :undoc-members:

    .. automethod:: aospy.data_loader.DataLoader.__init__

.. autoclass:: aospy.data_loader.DictDataLoader
    :members:
    :undoc-members:

    .. automethod:: aospy.data_loader.DictDataLoader.__init__

.. autoclass:: aospy.data_loader.NestedDictDataLoader
    :members:
    :undoc-members:

    .. automethod:: aospy.data_loader.NestedDictDataLoader.__init__

.. autoclass:: aospy.data_loader.GFDLDataLoader
    :members:
    :undoc-members:

    .. automethod:: aospy.data_loader.GFDLDataLoader.__init__

Variables and Regions
=====================

The ``Var`` and ``Region`` classes are used to represent,
respectively, physical quantities the user wishes to be able to
compute and geographical regions over which the user wishes to
aggregate their calculations.

Whereas the ``Proj`` - ``Model`` - ``Run`` hierarchy is used to
describe the data resulting from particular model simulations, ``Var``
and ``Region`` represent the properties of generic physical entities
that do not depend on the underlying data.

Var
---

.. autoclass:: aospy.var.Var
    :members:
    :undoc-members:

    .. automethod:: aospy.var.Var.__init__

Region
------

.. autoclass:: aospy.region.Region
    :members:
    :undoc-members:

    .. automethod:: aospy.region.Region.__init__

Calculations
============

``Calc`` is the engine that combines the user's specifications of (1)
the data on disk via ``Proj``, ``Model``, and ``Run``, (2) the
physical quantity to compute and regions to aggregate over via ``Var``
and ``Region``, and (3) the desired date range, time reduction method,
and other characteristics to actually perform the calculation

Whereas ``Proj``, ``Model``, ``Run``, ``Var``, and ``Region`` are all
intended to be saved in ``.py`` files for reuse, ``Calc`` objects are
intended to be generated dynamically by a main script and then not
retained after they have written their outputs to disk following the
user's specifications.

Moreover, if the ``main.py`` script is used to execute calculations,
no direct interfacing with ``Calc`` or it's helper class,
``CalcInterface`` is required by the user, in which case this section
should be skipped entirely.

Also included is the ``automate`` module, which enables aospy e.g. in
the main script to find objects in the user's object library that the
user specifies via their string names rather than having to import the
objects themselves.

CalcInterface and Calc
----------------------

.. autoclass:: aospy.calc.CalcInterface
    :members:
    :undoc-members:

    .. automethod:: aospy.calc.CalcInterface.__init__

.. autoclass:: aospy.calc.Calc
    :members:
    :undoc-members:

    .. automethod:: aospy.calc.Calc.__init__

automate
--------

.. automodule:: aospy.automate
    :members:
    :undoc-members:

operator
--------

.. warning::

    The ``operator`` module is in the process of being re-vamped and
    is therefore currently not supported.

.. automodule:: aospy.operator
    :members:
    :undoc-members:

Units and Constants
===================

aospy provides the classes ``Constant`` and ``Units`` for
representing, respectively, physical constants (e.g. Earth's
gravitational acceleration at the surface = 9.81 m/s^2) and physical
units (e.g. meters per second squared in that example).

aospy comes with several commonly used constants saved within the
``constants`` module in which the ``Constant`` class is also defined.
In contrast, there are no pre-defined ``Units`` objects; the user must
define any ``Units`` objects they wish to use (e.g. to populate the
``units`` attribute of their ``Var`` objects).

Similarly, whereas these baked-in ``Constant`` objects are used by
aospy in various places, aospy currently does not actually use the
``Var.units`` attribute or the ``Units`` class more generally; they
are for the user's own informational purposes.

constants
---------

.. automodule:: aospy.constants
    :members:
    :undoc-members:

units
-----

.. automodule:: aospy.units
    :members:
    :undoc-members:

.. note::

   There has been discussion of implementing units-handling upstream
   within xarray (see `here
   <https://github.com/pydata/xarray/issues/525>`_).  If and when that
   happens, the ``Units`` class will likely be deprecated and replaced
   with the upstream version.

Utilities
=========

aospy includes a number of utility functions that are used internally
and may also be useful to users for their own purposes.  These include
functions pertaining to input/output (IO), time arrays, andvertical
coordinates.

utils.io
--------

.. automodule:: aospy.utils.io
    :members:
    :undoc-members:

utils.times
-----------

.. automodule:: aospy.utils.times
    :members:
    :undoc-members:

utils.vertcoord
---------------

.. automodule:: aospy.utils.vertcoord
    :members:
    :undoc-members:
