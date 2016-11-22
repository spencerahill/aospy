.. aospy documentation master file, created by
   sphinx-quickstart on Wed Oct 12 16:24:08 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

aospy: automated climate data analysis and management
=====================================================
aospy is a Python-based tool for automating computations involving gridded climate data and the management of the results of those computations.

Motivations
-----------
Climate models generally output a wide array of useful quantities, but almost invariably not all needed quantities are directly outputted.  Even for those that are, further slicing and dicing in time and/or space are required.  Moreover, these multiple computations and spatiotemporal reductions are needed for not just a single simulation but across multiple models, simulations, time durations, subsets of the annual cycle, and so on.  Performing these computations across all of the desired parameter combinations quickly becomes impractical without some automation.  But even once some automation is in place, the resulting plethora of data quickly becomes unusable unless it is easily found and imbued with sufficient metadata to describe precisely what it is and how it was computed.

What aospy does
---------------
aospy provides functionality that enables users to perform commonly needed tasks in climate, weather, and related sciences as:

* repeating a calculation across multiple simulations in a single climate model
* repeating a calculation across the same simulations performed in multiple models, even if variable names or other quantities differ among the models or simulations
* repeating a calculation across multiple timespans, both in terms of the start date and end date and in terms of sub-annual sampling: e.g. annual mean, seasonal-means, monthly means.
* Computing multiple statistical  (e.g. mean, standard deviation) and physical (e.g. column integrals or averages; zonal integrals or averages) reductions on any given computation, both on a gridpoint-by-gridpoint basis and over an arbitrary number of geographical regions
* any combination of the above, plus many more!

Design philosophy
-----------------
Key to enabling this automation and handling of model- and simulation-level idiosyncrasies is separating

1. Code that describes the data that you want to work with
2. Code that specifies any physical calculations you eventually want to perform
3. Code that specifies the set of computations the user wishes to perform at a given time

For (1), the user defines objects at three distinct levels: `Proj`, `Model`, and `Run`, that specify where the data is located that you want to work with.  For (2), the user defines `Var` objects that describe the physical quantities to be computed, including any functions that transform one or more directly model-outputted quantities into the ultimately desired quantity.  Once these objects have been defined, the user can proceed with (3) via a simple script that specifies any models, simulations, physical quantities, etc. to be performed.

The run script can be modified and re-submitted as further calculations are desired.  Similarly, new objects can be defined at any time describing new simulations, models, or variables.

How to use
==========

Documentation
-------------
Formal documentation is available on ReadTheDocs: `<http://aospy.readthedocs.io/en/latest/>`_

*Note: A major improvement in the documentation (as well as a major refactor of aospy's internals) is currently under way.*

*Also note: aospy is in the alpha stage of development and therefore the code base is pretty unstable, both the internals and user-facing components.  We do not guarantee backwards compatibility of any features for future releases.*

Installation
------------
Simply clone the main Github repo::

  git clone https://www.github.com/spencerahill/aospy/aospy.git
  cd aospy
  python setup.py install

Once we reach the v0.1 release, aospy will become available on PyPI (i.e. ``pip install aospy``) and conda (``conda install aospy``).  For now, just use the above ``git clone`` method.  Please open an issue if you have any installation problems.

Examples
--------
The best way to understand the aospy workflow is to look at the `example IPython notebook <https://github.com/spencerahill/aospy/blob/develop/examples/example.ipynb>`_, which walks you through a simple use-case using example data that is also provided with the package.

Troubleshooting
---------------
Questions of any kind are welcome and best placed as Issues on the Github repo.

Development
===========
aospy is currently developed by `Spencer Hill <https://github.com/spencerahill>`_ and `Spencer Clark <https://github.com/spencerkclark>`_.  We are actively seeking new users and/or contributors, including during the ongoing refactor/documentation sprint.  Please get in touch!

API Reference
=============

.. toctree::
   :maxdepth: 2

   aospy

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
