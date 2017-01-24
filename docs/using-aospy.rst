###########
Using aospy
###########

This section provides a high-level summary of how to use aospy.  See
the Examples section and associated Jupyter Notebook for concrete
examples.

.. note::

   There is a non-trivial amount of effort required (mainly in
   creating and populating your object library, described below)
   before you will be able to perform any calculations.  However, once
   the object library is in place, there is essentially no limit to
   the number of calculations that can be performed on your data
   either together at the same time or at different times.  In other
   words, the spinup time should be well worth it.

Your aospy object library
=========================

The first step is writing the code that describes your data and the
quantities you eventually want to compute using it.  We refer to this
code collectively as your "object library".

Describing your data on disk
----------------------------

aospy needs to know where the data you want to use is located on disk
and how it is organized across different simulations, models, and
projects.  This involves a hierarchy of three classes, ``Proj``,
``Model``, and ``Run``.

1. ``Proj``: This represents a single project that involves analysis of
   data from one or more models and simulations.

2. ``Model``: This represents a single climate model, other numerical
   model, observational data source, etc.

3. ``Run``: This represents a single simulation, version of
   observational data, etc.

So each user's object library will contain one or more ``Proj``
objects, each of which will have one or more child ``Model`` objects,
which in turn will each have one or more child ``Run`` objects.

.. note::

   Currently, the Proj-Model-Run hierarchy is rigid, in that each Run
   has a parent Model, and each Model has a parent Proj.  For small
   projects, this can lead to a lot of boilerplate code.  Work is
   ongoing to relax this constraint to facilitate easier exploratory
   analysis.

Physical variables
------------------

The ``Var`` class is used to represent physical variables,
e.g. precipitation rate or potential temperature.  This includes both
variables which are directly available in netCDF files (e.g. they were
directly outputted by your climate model) as well as those fields that
must be computed from other variables (e.g. they weren't directly
outputted but can be computed from other variables that were
outputted).

Geographical regions
--------------------

The ``Region`` class is used to define geographical regions over which
quantities can be averaged (in addition to gridpoint-by-gridpoint
values).  Like ``Var`` objects, they are more generic than the objects
of the ``Proj`` - ``Model`` - ``Run`` hierarchy, in that they
correspond to the generic physical quantities/regions rather than the
data of a particular project, model, or simulation.

Configuring your object library
===============================

Required components
-------------------

In order for your object library to work with the main script, it must
include the following two objects:

1. ``projs`` : A container of ``Proj`` objects
2. ``variables`` : A container of ``Var`` objects

(The ``Model``, ``Run``, and ``Region`` objects are all included
within their parent ``Proj`` objects and thus don't require analogous
top-level containers.)

These must be accessible from the object library's toplevel namespace,
i.e. the Python commands ``import my_obj_lib.projs`` and ``import
my_obj_lib.variables`` must work, where ``my_obj_lib`` is the name
you've given to your library.  Which leads to the next topic: how to
structure your object library within one or more ``.py`` files.

File/directory structure
------------------------

The simplest way to structure your object library is to define
everything in a single module (i.e. a single ``.py`` file).  This
works great for small projects and for initially trying out aospy.

As an object library grows, however, it can become desirable to split
it into multiple ``.py`` files.  This effectively changes it from a
module to a proper Python package.  Python packages require a specific
directory structure and specification of things to include at each
level via ``__init__.py`` files.  See the `official documentation
<https://docs.python.org/3.6/tutorial/modules.html#packages>`_ on
packages for further guidance.

For an example of a large object library that is structured as a
proper package, see `here
<https://github.com/spencerahill/aospy-obj-lib>`_.

Making your object library visible to Python
--------------------------------------------

Whether it is structured as a single module or as a proper package,
you'll likely have to add the directory containing your object library
to the ``PYTHONPATH`` environment variable in order for Python to be
able to import it::

  export PYTHONPATH=/path/to/your/object/library:${PYTHONPATH}

Of course, replace ``/path/to/your/object/library`` with the actual
path to yours.  This command places your object library at the front
of the ``PYTHONPATH`` environment variable, which is essentially the
first place where Python looks to find packages and modules to be
imported.

.. note::

   It's convenient to copy this command into your shell profile (e.g.,
   for the bash shell on Linux or Mac, ``~/.bash_profile``) so that
   you don't have to call it again in every new terminal session.

.. note::

   For object libraries structured as packages, it is also possible to
   properly install your object library by creating a properly set-up
   ``setup.py`` file and ``python setup.py install``.  But unless
   you're prevented from modifying ``PYTHONPATH`` for some reason,
   there's no advantage of this versus the simpler
   ``PYTHONPATH`` alternative above.

Once this has been done, you should be able to import your object
library from within Python via ``import my_obj_lib``, where
``my_obj_lib`` is the name you've given to your library.  You will not
be able to use the main script until this works.

Executing calculations
======================

The main script contents
------------------------

Calculations are performed by specifying in a "main script" the
desired parameters and then running the script.

We provide a template main script within aospy.  You should copy it to
the location of your choice and in the copy replace the given names
with the names of your own project, model, etc. objects that you want
to perform computations on.  (If you accidentally change the original,
you can always get a `fresh copy from Github
<https://github.com/spencerahill/aospy/tree/develop/examples>`_.)

Except where noted otherwise in the template script's comments, all
parameters should be submitted as lists, even if they are a single
element.  E.g. ``models = ['name-of-my-model']``.

.. note::

   Although the main script is the recommended way to perform
   calculations, it's possible to submit calculations by other means.
   For example, one could explicitly create ``Calc`` objects and call
   their ``compute`` method, as is done in the example Jupyter
   notebook.

Running the main script
-----------------------
Once the main script parameters are all modified as desired, execute
the script from the command line as follows ::

  /path/to/your/main.py

This should generate a text summary of the specified parameters and a
prompt as to whether to proceed or not with the calculations.  An
affirmative response then triggers the calculations to execute.

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

As the calculations are performed, logging information will be printed
to the terminal displaying their progress.

Parallelized calculations
-------------------------

The calculations generated by the main script can be executed in
parallel provided the optional dependency ``multiprocess`` is
installed.  (It is available via pip: ``pip install multiprocess``.)
Otherwise, or if the user sets ``parallelize`` to ``False`` in the main
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

See the API reference documentation of ``CalcInterface`` for explanation of each of these components of the path and file name.

Under the hood
--------------

The main script encodes each permutation of the input parameters into
a ``CalcInterface`` object.  This object, in turn, is used to
instantiate a ``Calc`` object.  The ``Calc`` object, in turn, performs
the calculation.

Unlike ``Proj``, ``Model``, ``Run``, ``Var``, and ``Region``, these
objects are not intended to be saved in ``.py`` files for continual
re-use.  Instead, they are generated as needed, perform their desired
tasks, and then go away.

See the API reference documentation for further details.
