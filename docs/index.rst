.. aospy documentation master file, created by
   sphinx-quickstart on Wed Oct 12 16:24:08 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: aospy_logo.png
   :alt: aospy logo
   :align: center
   :height: 300px
   :width: 300px
   :name: aospy-logo

#####################################################
aospy: automated climate data analysis and management
#####################################################

**aospy** is an open source Python package for automating computations
that use gridded climate data (namely data stored as netCDF files)
and the management of the results of those computations.

Once a user describes where their data is stored on disk using aospy's
built-in tools, they can subsequently use the provided main script at
any time to fire off calculations to be performed in parallel using
the permutation of an arbitrary number of climate models, simulations,
variables to be computed, date ranges, sub-annual-sampling, and many
other parameters.  Their results get saved in a highly structured
directory format as netCDF files.

The eventual goal is for aospy to become the "industry standard" for
gridded climate data analysis and, in so doing, accelerate progress in
climate science and make the results of climate research more easily
reproducible and shareable.  aospy relies heavily on the `xarray
<http://xarray.pydata.org>`_ package.

Documentation
=============
.. toctree::
   :maxdepth: 1

   whats-new
   overview
   using-aospy
   examples
   install
   api

Get in touch
============

- Troubleshooting: We are actively seeking new users and are eager to
  help you get started with aospy!  Usage questions, bug reports, and
  any other correspondence are all welcome and best placed as `Issues
  on the Github repo <https://github.com/spencerahill/aospy>`_.
- Contributing: We are also actively seeking new developers!  Please
  get in touch by opening an Issue or submitting a Pull Request.

License
=======

aospy is freely available under the open source `Apache License
<http://www.apache.org/licenses/>`_.

History
=======

aospy was originally created by Spencer Hill as a means of automating
calculations required for his Ph.D. work across many climate models
and simulations.  Starting in 2014, Spencer Clark became aospy's
second user and developer.  The first official release on PyPI was
v0.1 on January 24, 2017.
