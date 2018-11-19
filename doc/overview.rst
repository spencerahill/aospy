.. _overview:

####################
Overview: Why aospy?
####################

Use cases
=========

If you've ever found yourself saying or thinking anything along these
lines, then aospy may be a great tool for you:

- "What is this precip01.dat file in my home directory that I vaguely
  remember creating a few weeks ago?  Better just delete it and re-do
  the calculation to be sure."
- "I really need to analyze variable X.  But the models/observational
  products I'm using don't provide it.  They do provide variables Y
  and Z, from which I can compute quantity X."
- "I need to calculate quantity X from 4 different simulations
  repeated in 5 different models; computing monthly, seasonal, and
  annual averages and standard deviations; gridpoint-by-gridpoint and
  averaged over these 10 regions.  That's impractical -- I'll just do
  a small subset."

Each of these common problems is easily solved by using aospy.

.. _design-philosophy:

Design philosophy
=================

aospy's ability to automate calculations (while properly handling
model- and simulation-level idiosyncrasies) relies on separating your
code into three distinct categories.

1. Code characterizing the data you want to work with: "Where is your
   data and what does it represent?"
2. Code describing abstract physical quantities and geographical
   regions you eventually want to examine: "What things are you
   generally interested in calculating?"
3. Code specifying the exact parameters of calculations you want to
   perform right now: "Ok, time to actually crunch some numbers.
   Exactly what all do you want to compute from your data, and how do
   you want to slice and dice it?"

How you'll actually interact with aospy in order to achieve each of
these three steps is described in the :ref:`using-aospy` section of
this documentation, and explicit examples using included sample data
are in the :ref:`examples` section.

Open Science & Reproducible Research
====================================

This separation of your code into three categories promotes `open
science <https://en.wikipedia.org/wiki/Open_science>`_ and
`reproducible research
<https://en.wikipedia.org/wiki/Reproducibility#Reproducible_research>`_
in multiple ways:

- By separating code that describes where your data is from the code
  used to compute particular physical quantities, the latter can be
  written in a generic form that closely mimics the physical form of
  the particular expression.  The resulting code is easier to read and
  debug and therefore to share with others.
- By enabling automation of calculations across an arbitrary number of
  parameter combinations, aospy facilitates more rigorous analyses
  and/or analyses encompassing a larger span of input data than would
  otherwise be possible.
- By outputting the results of calculations as netCDF files in a
  highly organized directory structure with lots of metadata embued
  within the file path, file name, and the netCDF file itself, aospy
  facilitates the sharing of results with others (including your
  future self that has forgotten the myriad details of how you
  have computed things right now).
