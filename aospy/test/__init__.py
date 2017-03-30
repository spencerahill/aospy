import pytest

try:
    import IPython
    has_ipython = True
except ImportError:
    has_ipython = False

try:
    import runipy
    has_runipy = True
except ImportError:
    has_runipy = False

try:
    import matplotlib
    has_matplotlib = True
except ImportError:
    has_matplotlib = False


def requires_ipython(test):
    return test if has_ipython else pytest.mark.skip('requires IPython')(test)

def requires_runipy(test):
    return test if has_runipy else pytest.mark.skip('requires runipy')(test)

def requires_matplotlib(test):
    if has_matplotlib:
        return test
    else:
        return pytest.mark.skip('requires matplotlib')(test)
