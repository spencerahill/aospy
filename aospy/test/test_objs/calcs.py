from aospy.utils import dp_from_ps
from aospy.__config__ import PFULL_STR


def dp(ps, bk, pk, arr):
    """Pressure thickness of hybrid coordinate levels from surface pressure."""
    return dp_from_ps(bk, pk, ps, arr[PFULL_STR])
