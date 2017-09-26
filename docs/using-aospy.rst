.. _using-aospy:

###########
Using aospy
###########

This section provides a high-level summary of how to use aospy.  See
the :ref:`Overview <overview>` section of this documentation for more
background information, or the :ref:`Examples` section for concrete
examples.

Your aospy object library
=========================

The first step is writing the code that describes your data and the
quantities you eventually want to compute using it.  We refer to this
code collectively as your "object library".

Describing your data on disk
----------------------------

aospy needs to know where the data you want to use is located on disk
and how it is organized across different projects, models, and model
runs (i.e. simulations).  This involves a hierarchy of three classes,
:py:class:`aospy.Proj`, :py:class:`aospy.Model`, and
:py:class:`aospy.Run`.

1. :py:class:`aospy.Proj`: This represents a single project that
   involves analysis of data from one or more models and simulations.

2. :py:class:`aospy.Model`: This represents a single climate model,
   other numerical model, observational data source, etc.

3. :py:class:`aospy.Run`: This represents a single simulation,
   version of observational data, etc.

So each user's object library will contain one or more
:py:class:`aospy.Proj` objects, each of which will have one or more
child :py:class:`aospy.Model` objects, which in turn will each have
one or more child :py:class:`aospy.Run` objects.

.. note::

   Currently, the Proj-Model-Run hierarchy is rigid, in that each Run
   has a parent Model, and each Model has a parent Proj.  Work is
   ongoing to relax this to a more generic parent-child framework.

Physical variables
------------------

The :py:class:`aospy.Var` class is used to represent physical variables,
e.g. precipitation rate or potential temperature.  This includes both
variables which are directly available in netCDF files (e.g. they were
directly outputted by your model or gridded data product) as well as
those fields that must be computed from other variables (e.g. they
weren't directly outputted but can be computed from other variables
that were outputted).

Geographical regions
--------------------

The :py:class:`aospy.Region` class is used to define geographical
regions over which quantities can be averaged (in addition to
gridpoint-by-gridpoint values).  Like :py:class:`aospy.Var` objects,
they are more generic than the objects of the :py:class:`aospy.Proj` -
:py:class:`aospy.Model` - :py:class:`aospy.Run` hierarchy, in that
they correspond to the generic physical quantities/regions rather than
the data of a particular project, model, or simulation.

Object library structure
------------------------

The officially supported way to submit calculations is the
:py:meth:`aospy.submit_mult_calcs` function.  In order for this to
work, your object library must follow one or the other of these
structures:

1. All :py:class:`aospy.Proj` and :py:class:`aospy.Var` objects are
   accessible as attributes of your library.  This means that
   ``my_obj_lib.my_obj`` works, where ``my_obj_lib`` is
   your object library, and ``my_obj`` is the object in question.
2. All :py:class:`aospy.Proj` objects are stored in a container called
   ``projs``, where ``projs`` is an attribute of your library
   (i.e. ``my_obj_lib.projs``).  And likewise for
   :py:class:`aospy.Var` objects in a ``variables`` attribute.

Beyond that, you can structure your object library however you wish.
In particular, it can be structured as a Python module (i.e. a single
".py" file) or as a package (i.e. multiple ".py" files linked
together; see the `official documentation
<https://docs.python.org/3.6/tutorial/modules.html#packages>`_ on
package structuring).

A single module works great for small projects and for initially
trying out aospy (this is how the example object library,
:py:mod:`aospy.examples.example_obj_lib`, is structured).  But as
your object library grows, it can become easier to manage as a package
of multiple files.  For an example of a large object library that is
structured as a formal package, see `here
<https://github.com/spencerahill/aospy-obj-lib>`_.

Accessing your library
----------------------

If your current working directory is the one containing your library,
you can import your library via ``import my_obj_lib`` (replacing
``my_obj_lib`` with whatever you've named yours) in order to pass it
to :py:meth:`aospy.submit_mult_calcs`.

Once you start using aospy a lot, however, this requirement of being
in the same directory becomes cumbersome.  As a solution, you can add
the directory containing your object library to the ``PYTHONPATH``
environment variable.  E.g if you're using the bash shell: ::

  export PYTHONPATH=/path/to/your/object/library:${PYTHONPATH}

Of course, replace ``/path/to/your/object/library`` with the actual
path to yours.  This command places your object library at the front
of the ``PYTHONPATH`` environment variable, which is essentially the
first place where Python looks to find packages and modules to be
imported.  (For more, see Python's `official documentation on
PYTHONPATH <https://docs.python.org/3.6/using/cmdline.html>`_).

.. note::

   It's convenient to copy this command into your shell profile (e.g.,
   for the bash shell on Linux or Mac, ``~/.bash_profile``) so that
   you don't have to call it again in every new terminal session.

To test this is working, run ``python -c "import my_obj_lib"`` from a
directory other than where the library is located (again replacing
``my_obj_lib`` with the name you've given to your library).  If this
runs without error, you should be good to go.

Executing calculations
======================

As noted above, the officially supported way to submit calculations is the
:py:meth:`aospy.submit_mult_calcs` function.

We provide a template "main" script with aospy that uses this
function.  We recommend copying it to the location of your choice.  In
the copy, replace the example object library and associated objects
with your own.  (If you accidentally change the original, you can
always get a `fresh copy from Github
<https://github.com/spencerahill/aospy/blob/develop/aospy/examples/aospy_main.py>`_).

Running the main script
-----------------------
Once the main script parameters are all modified as desired, execute
the script from the command line as follows ::

  /path/to/your/aospy_main.py

This should generate a text summary of the specified parameters and a
prompt as to whether to proceed or not with the calculations.  An
affirmative response then triggers the calculations to execute.

.. note::

   You may need to change the permissions on the file to make it
   executable.  E.g. from a Mac or Linux: `chmod u+x
   /path/to/your/aospy_main.py`.  Alternatively you can call python or
   IPython from the command line to run it: `python
   /path/to/your/aospy_main.py` or `ipython /path/to/your/aospy_main.py`.

Specifically, the parameters are permuted over all possible
combinations.  So, for example, if two model names and three variable
names were listed and all other parameters had only one element, six
calculations would be generated and executed.  There is no limit to
the number of permutations.

.. note::

   You can also call the main script interactively within an IPython
   session via ``%run /path/to/your/main.py`` or, from the command
   line, run the script and then start an interactive IPython session
   via ``ipython -i /path/to/your/main.py``.

   Or you can call :py:func:`aospy.submit_mult_calcs` directly within
   an interactive session.

As the calculations are performed, logging information will be printed
to the terminal displaying their progress.

Parallelized calculations
-------------------------

The calculations generated by the main script can be executed in
parallel using ``dask.distributed``. aospy will either automatically
set up a ``dask.distributed.LocalCluster`` to perform the calculations,
or one can optionally specify an external ``distributed.Client`` to delegate
the work. Otherwise, or if the user sets ``parallelize=False`` in the
``calc_exec_options`` argument of :py:func:`aospy.submit_mult_calcs`,
script, the calculations will be executed one-by-one.

Particularly on instititutional clusters with many cores, this
parallelization yields an impressive speed-up when multiple
calculations are generated.

.. note::

   When calculations are performed in parallel, often the logging
   information from different calculations running simultameously end
   up interwoven with one another, leading to output that is confusing
   to follow.  Work is ongoing to improve the logging output when the
   computations are parallelized.

Finding the output
------------------

aospy saves the results of all calculations as netCDF files and embeds
metadata describing it within the netCDF files, in their filenames,
and in the directory structure within which they are saved.

- Directory structure:
  ``/path/to/aospy-rootdir/projname/modelname/runname/varname``
- File name :
  ``varname.intvl_out.dtype_out_time.'from_'intvl_in'_'dtype_in_time.model.run.date_range.nc``

See the :ref:`api-ref` on :py:class:`aospy.CalcInterface` for
explanation of each of these components of the path and file name.

Under the hood
==============

:py:func:`aospy.submit_mult_calcs` creates a :py:class:`aospy.CalcSuite`
object that permutes over the provided lists of calculation
specifications, encoding each permutation into a
:py:class:`aospy.CalcInterface` object.

.. note::

   Actually, when multiple regions and/or output time/regional
   reductions are specified, these all get passed to each
   :py:class:`aospy.CalcInterface` object rather than being permuted
   over.  They are then looped over during the subsequent
   calculations.  This is to prevent unnecessary re-loading and
   re-computing, because, for a given simulation/variable/etc., all
   regions and reduction methods use the same data.

Each :py:class:`aospy.CalcInterface` object, in turn, is used to
instantiate a :py:class:`aospy.Calc` object.  The
:py:class:`aospy.Calc` object, in turn:

- loads the required netCDF data given its simulation, variable, and date range
- (if necessary) further truncates the data in time (i.e. to the given
  subset of the annual cycle, and/or if the requested date range
  doesn't exactly align with the time chunking of the input netCDF
  files)
- (if the variable is a function of other variables) executes the
  function that computes the calculation using this loaded and
  truncated data
- applies all specified temporal and regional time reductions
- writes the results (plus additional metadata) to disk as netCDF
  files and appends it to its own ``data_out`` attribute

.. note::

   Unlike :py:class:`aospy.Proj`, :py:class:`aospy.Model`,
   :py:class:`aospy.Run`, :py:class:`aospy.Var`, and
   :py:class:`aospy.Region`, these objects are not intended to be
   saved in ``.py`` files for continual re-use.  Instead, they are
   generated as needed, perform their desired tasks, and then go away.

See the :ref:`API reference <api-ref>` documentation for further details.
