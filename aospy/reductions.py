"""Functionality relating to spatiotemporal reduction operations."""
from .internal_names import TIME_STR, TIME_WEIGHTS_STR, YEAR_STR
from .utils.times import yearly_average


class Reduction(object):
    """Base class for spatiotemporal reductions."""
    instances = {}

    @classmethod
    def register_reduction(cls, reduction):
        label = reduction.label
        if label in cls.instances:
            raise KeyError("A Reduction object already exists with the "
                           "desired label '{0}', and duplicates are not "
                           "allowed. The existing object is: "
                           "{1}".format(label, cls.instances[label]))
        else:
            cls.instances[label] = reduction

    def __init__(self, reduction_func, label):
        self.register_reduction(self)
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


def _in_year_avg(arr):
    return yearly_average(arr, arr[TIME_WEIGHTS_STR])


no_reduc = Reduction(lambda x: x, 'full')
time_avg = Reduction(lambda x: x.mean(TIME_STR), 'av')
yearly_ts = Reduction(lambda x: _in_year_avg(x), 'ts')
yearly_stdev = Reduction(lambda x: _in_year_avg(x).std(YEAR_STR), 'std')
