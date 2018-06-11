"""Functionality relating to spatiotemporal reduction operations."""
from .internal_names import YEAR_STR


_TIME_DEFINED_REDUCTIONS = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']


class Reduction(object):
    """Base class for spatiotemporal reductions."""
    def __init__(self, reduction_func, label):
        self._reducer = reduction_func
        if isinstance(label, str):
            self._label = label
        else:
            raise TypeError("'label' must be a string; got type "
                            "{}".format(label))

    @property
    def reduce(self):
        return self._reducer

    @reduce.setter
    def reduce(self, value):
        raise ValueError("'reduce' method cannot be modified after Reduction "
                         "object has been created.")

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        raise ValueError("'label' method cannot be modified after Reduction "
                         "object has been created.")

    def __call__(self, data):
        return self.reduce(data)


no_reduc = Reduction(lambda x: x)
time_avg = Reduction(lambda x: x.mean(YEAR_STR))
time_stdev = Reduction(lambda x: x.std(YEAR_STR))
