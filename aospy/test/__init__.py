import pytest


def requires_pytest_catchlog(test):
    try:
        import pytest_catchlog
    except ImportError:
        return pytest.mark.skip('requires pytest_catchlog')(test)
    else:
        return test
