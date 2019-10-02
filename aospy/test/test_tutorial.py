import os
import sys
import warnings

import pytest

import aospy


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
