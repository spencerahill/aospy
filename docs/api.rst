.. _api-ref:

#############
API Reference
#############

Here we provide the reference documentation for aospy's public API.
If you are new to the package and/or just trying to get a feel for the
overall workflow, you are better off starting in the :ref:`Overview
<overview>`, :ref:`using-aospy`, or :ref:`examples` sections of this
documentation.

.. warning::

   aospy is under active development.  While we strive to maintain
   backwards compatibility, it is likely that some breaking changes to
   the codebase will occur in the future as aospy is improved.

Core Hierarchy for Input Data
=============================

aospy provides three classes for specifying the location and
characteristics of data saved on disk as netCDF files that the user
wishes to use as input data for aospy calculations: :py:class:`Proj`,
:py:class:`Model`, and :py:class:`Run`.

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

:py:class:`Run` objects rely on a helper "data loader" to specify how to find
their underlying data that is saved on disk.  This mapping of
variables, time ranges, and potentially other parameters to the
location of the corresponding data on disk can differ among modeling
centers or even between different models at the same center.

Currently supported data loader types are :py:class:`DictDataLoader`,
:py:class:`NestedDictDataLoader`, and :py:class:`GFDLDataLoader` Each of these inherit
from the abstract base :py:class:`DataLoader` class.

.. note::

   An important consideration can be the datatype used to store values in
   your datasets.  In particular, if the float32 datatype is used in
   storage, it can lead to undesired inaccuracies in the computation of
   reduction operations (like means) due to upstream issues (see
   `pydata/xarray#1346 <https://github.com/pydata/xarray/issues/1346>`_ for
   more information).  To address this it is recommended to always upcast
   float32 data to float64.  This behavior is turned on by default.  If you
   would like to disable this behavior you can set the ``upcast_float32``
   argument in your ``DataLoader`` constructors to ``False``.

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

The :py:class:`Var` and :py:class:`Region` classes are used to represent,
respectively, physical quantities the user wishes to be able to
compute and geographical regions over which the user wishes to
aggregate their calculations.

Whereas the :py:class:`Proj` - :py:class:`Model` - :py:class:`Run` hierarchy is used to
describe the data resulting from particular model simulations, :py:class:`Var`
and :py:class:`Region` represent the properties of generic physical entities
that do not depend on the underlying data.

Var
---

.. autoclass:: aospy.var.Var
    :members:
    :undoc-members:

    .. automethod:: aospy.var.Var.__init__

.. note::
    While for the sake of tracking metadata we encourage users to add a ``units``
    attribute to each of their :py:class:`Var` objects, these ``units`` attributes provide
    nothing more than descriptive value.  One day we hope DataArrays
    produced by loading or computing variables will be truly units-aware
    (e.g. adding or subtracting two DataArrays with different units will lead to
    an error, or multiplying two DataArrays will result in a DataArray with new units),
    but we will leave that to upstream libraries to implement
    (see `pydata/xarray#525 <https://github.com/pydata/xarray/issues/525>`_
    for more discussion).

Region
------

.. autoclass:: aospy.region.Region
    :members:
    :undoc-members:

    .. automethod:: aospy.region.Region.__init__

Calculations
============

:py:class:`Calc` is the engine that combines the user's specifications of (1)
the data on disk via :py:class:`Proj`, :py:class:`Model`, and :py:class:`Run`, (2) the
physical quantity to compute and regions to aggregate over via :py:class:`Var`
and :py:class:`Region`, and (3) the desired date range, time reduction method,
and other characteristics to actually perform the calculation

Whereas :py:class:`Proj`, :py:class:`Model`, :py:class:`Run`, :py:class:`Var`, and :py:class:`Region` are all
intended to be saved in ``.py`` files for reuse, :py:class:`Calc` objects are
intended to be generated dynamically by a main script and then not
retained after they have written their outputs to disk following the
user's specifications.

Moreover, if the ``main.py`` script is used to execute calculations,
no direct interfacing with :py:class:`Calc` or it's helper class,
:py:class:`CalcInterface` is required by the user, in which case this section
should be skipped entirely.

Also included is the :py:class:`automate` module, which enables aospy e.g. in
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
