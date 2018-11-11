import os
import sys
import warnings

import pytest

import aospy


ON_WINDOWS = sys.platform == 'win32' or sys.platform == 'win64'


@pytest.mark.skipif(
    ON_WINDOWS, reason='This method of running Jupyter '
    'notebooks on AppVeyor no longer seems to work.')
def test_tutorial_notebook():
    pytest.importorskip('nbformat')
    pytest.importorskip('nbconvert')
    pytest.importorskip('matplotlib')

    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor

    rootdir = os.path.join(aospy.__path__[0], 'examples')
    with open(os.path.join(rootdir, 'tutorial.ipynb')) as nb_file:
        notebook = nbformat.read(nb_file, as_version=nbformat.NO_CONVERT)
    kernel_name = 'python' + str(sys.version[0])
    ep = ExecutePreprocessor(kernel_name=kernel_name)
    with warnings.catch_warnings(record=True):
        ep.preprocess(notebook, {})
