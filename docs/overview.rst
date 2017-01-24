Overview: Why aospy?
====================

Motivations
-----------

Climate models generally output a wide array of useful quantities, but
almost invariably not all needed quantities are directly outputted.
Even for those that are, further slicing and dicing in time and/or
space are required.  Moreover, these multiple computations and
spatiotemporal reductions are needed for not just a single simulation
but across multiple models, simulations, time durations, subsets of
the annual cycle, and so on.

Performing these computations across all of the desired parameter
combinations quickly becomes impractical without some automation.  But
even once some automation is in place, the resulting plethora of data
quickly becomes unusable unless it is easily found and imbued with
sufficient metadata to describe precisely what it is and how it was
computed.

What aospy does
---------------

aospy provides functionality that enables users to perform commonly
needed tasks in climate, weather, and related sciences as:

* Repeating a calculation across multiple simulations in a single
  climate model
* Repeating a calculation across the same simulations performed in
  multiple models, even if variable names or other quantities differ
  among the models or simulations
* Repeating a calculation across multiple timespans, both in terms of
  the start date and end date and in terms of sub-annual sampling:
  e.g. annual mean, seasonal-means, monthly means.
* Computing multiple statistical (e.g. mean, standard deviation) and
  physical (e.g. column integrals or averages; zonal integrals or
  averages) reductions on any given computation, both on a
  gridpoint-by-gridpoint basis and over an arbitrary number of
  geographical regions
* Any combination of the above, plus many more!

Design philosophy
-----------------

Key to enabling this automation and handling of model- and
simulation-level idiosyncrasies is separating

1. Code that describes the data that you want to work with
2. Code that specifies any physical calculations you eventually want
   to perform
3. Code that specifies the set of computations the user wishes to
   perform at a given time

For (1), the user defines objects at three distinct levels: `Proj`,
`Model`, and `Run`, that specify where the data is located that you
want to work with.  For (2), the user defines `Var` objects that
describe the physical quantities to be computed, including any
functions that transform one or more directly model-outputted
quantities into the ultimately desired quantity, as well as `Region`
objects that describe any geographical regions over which to perform
averages.  Once these objects have been defined, the user can proceed
with (3) via a simple script that specifies any models, simulations,
physical quantities, etc. to be performed.

The run script can be modified and re-submitted as further
calculations are desired.  Similarly, new objects can be defined at
any time describing new simulations, models, or variables.  More
detailed instructions are available in the "Using aospy" section.

Open Science & Reproducible Research
------------------------------------

aospy promotes `open science
<https://en.wikipedia.org/wiki/Open_science>`_ and `reproducible
research
<https://en.wikipedia.org/wiki/Reproducibility#Reproducible_research>`_
in multiple ways:

- By separating code that describes where your data is from the code
  used to compute particular physical quantities, the latter can be
  written in a generic form that closely mimics the physical form of
  the particular expression.  The resulting code is easier to read and
  debut and therefore to share with others.
- By enabling automation of calculations across an arbitrary number of
  parameter combinations, aospy facilitates more rigorous analyses
  and/or analyses encompassing a larger span of input data than would
  otherwise be possible.
- By outputting the results of calculations as netCDF files in a
  highly organized directory structure with lots of metadata embued
  within the file path, file name, and the netCDF file itself, aospy
  facilitates the sharing of data with others.

  It also enhances the usability of one's own data, providing a remedy
  to the familiar refrain among scientists along the lines of "What is
  this data1.txt file that was created six months ago?  Better just
  delete it and re-do the calculations to be sure.")
