Installation
============

Supported platforms
-------------------

- Operating system: Windows, MacOS/Mac OS X, and Linux
- Python: 2.7, 3.4, 3.5, and 3.6

.. note::
   For users that are new to Python, we recommend installing the free `Anaconda distribution <https://www.continuum.io/downloads>`_ of Python, which works on Mac, Linux, and Windows and both on normal computers and institutional clusters and doesn't require root permissions.

Recommended installation method: pip
------------------------------------

aospy is available from the official `Python Packaging Index (PyPI) <https://pypi.io>`_ via ``pip``::

  pip install aospy

.. note:: We are currently working on adding aospy to `conda-forge <https://conda-forge.github.io/>`_ such that it will be installable via `conda <http://conda.pydata.org/docs/>`_ (i.e. ``conda install aospy -c conda-forge``).

Alternative method: clone from Github
-------------------------------------

You can also clone the `Github repo <https://github.com/spencerahill/aospy>`_ ::

  git clone https://www.github.com/spencerahill/aospy/aospy.git
  cd aospy
  python setup.py install

Verifying proper installation
-----------------------------

Once installed via either method, you can run aospy's suite of tests using `py.test <http://doc.pytest.org/>`_.  From the top-level directory of the aospy installation ::

  pip install pytest  # if you don't have it already
  py.test aospy

If this results in any error messages or test failures, something has gone wrong.

Troubleshooting
---------------

Please open an `Issue on Github <https://github.com/spencerahill/aospy/issues/new>`_ if you have problems with any of these installation methods.

