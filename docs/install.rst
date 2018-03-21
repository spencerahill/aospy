.. _install:

Installation
============

Supported platforms
-------------------

- Operating system: Windows, MacOS/Mac OS X, and Linux
- Python: 2.7, 3.5, and 3.6

.. note::

   We highly recommend installing and using the free `Anaconda
   <https://www.continuum.io/downloads>`_ (or Miniconda) distribution
   of Python, which works on Mac, Linux, and Windows, both on normal
   computers and institutional clusters and doesn't require root
   permissions.

   See e.g. `this blog post
   <https://medium.com/@rabernat/custom-conda-environments-for-data-science-on-hpc-clusters-32d58c63aa95#.hqyl6y38i>`_
   for instructions on how to install Miniconda on a large cluster
   without root permission.

Recommended installation method: conda
--------------------------------------

The recommended installation method is via `conda
<http://conda.pydata.org/docs/>`_ ::

  conda install -c conda-forge aospy

Alternative method #1: pip
--------------------------

aospy is available from the official `Python Packaging Index (PyPI)
<https://pypi.io>`_ via ``pip``::

  pip install aospy

Alternative method #2: clone from Github
----------------------------------------

You can also directly clone the `Github repo
<https://github.com/spencerahill/aospy>`_ ::

  git clone https://www.github.com/spencerahill/aospy.git
  cd aospy
  python setup.py install

Verifying proper installation
-----------------------------

Once installed via any of these methods, you can run aospy's suite of
tests using `py.test <http://doc.pytest.org/>`_.  From the top-level
directory of the aospy installation ::

  conda install pytest  # if you don't have it already; or 'pip install pytest'
  py.test aospy

If you don't know the directory where aospy was installed, you can find it via ::

  python -c "import aospy; print(aospy.__path__[0])"

If the pytest command results in any error messages or test failures,
something has gone wrong, and please refer to the Troubleshooting
information below.

Troubleshooting
---------------

Please search through our `Issues page`_ on Github and our `mailing
list`_ to see if anybody else has had the same problem you're facing.
If none do, then please send a message on the list or open a new
Issue.

Please don't hesitate, especially if you're new to the package and/or
Python!  We are eager to help people get started.  Even you suspect
it's something simple or obvious -- that just means we didn't make it
clear/easy enough!

.. _Issues page: https://github.com/spencerahill/aospy/issues
.. _mailing list: https://groups.google.com/d/forum/aospy
