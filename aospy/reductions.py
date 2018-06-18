"""Functionality relating to spatiotemporal reduction operations."""
from .internal_names import TIME_STR, TIME_WEIGHTS_STR, YEAR_STR
from .utils.times import yearly_average


def _generic_register_reduction(cls, reduction, label):
    if label in cls.instances:
        raise KeyError("A Reduction object already exists with the "
                       "desired label '{0}', and duplicates are not "
                       "allowed. The existing object is: "
                       "{1}".format(label, cls.instances[label]))
    else:
        cls.instances[label] = reduction


class Reduction(object):
    """Base class for spatiotemporal reductions."""
    # TODO: make `instances` read-only except for from the
    # `_register_reduction` call within `__init__`
    instances = {}

    @classmethod
    def _register_reduction(cls, reduction, label):
        _generic_register_reduction(cls, reduction, label)

    def __init__(self, reduction_func, label):
        self._register_reduction(self, label)
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

    def __call__(self, data, *args, **kwargs):
        return self.reduce(data, *args, **kwargs)


class GriddedReduction(Reduction):
    # TODO: I think I am overwriting the parent class's `instances` attribute.
    # How can I register the instances of each subclass?  Or is this fine,
    # since the generic Reduction class won't ever be used, and other Reduction
    # subclasses won't overwrite each other's `instances` dictionary, just
    # the (empty) parent one.
    instances = {}

    @classmethod
    def _register_reduction(cls, reduction, label):
        _generic_register_reduction(cls, reduction, label)


class RegionReduction(Reduction):
    instances = {}

    @classmethod
    def _register_reduction(cls, reduction, label):
        _generic_register_reduction(cls, reduction, 'reg.' + label)

    def __init__(self, reduction_func, label, coord_reduction=None):
        self._register_reduction(self, label)
        self._reducer = reduction_func
        self._coord_reducer = coord_reduction

        if isinstance(label, str):
            self._label = label
        else:
            raise TypeError("'label' must be a string; got type "
                            "{}".format(label))

    @property
    def label(self):
        return 'reg.' + self._label

    @label.setter
    def label(self, value):
        raise ValueError("'label' method cannot be modified after Reduction "
                         "object has been created.")

    def __call__(self, data, region, *args, coord_name=None, coord_data=None,
                 **kwargs):
        data_out = self.reduce(region.ts(data))
        if self._coord_reducer is not None and coord_name is not None:
            if coord_data is None:
                coord_data = data[coord_name]
            data_out[coord_name] = self._coord_reducer(coord_data, *args,
                                                       **kwargs)
        return data_out


def _null_reduc(arr):
    return arr


def _time_avg(arr):
    return arr.mean(TIME_STR)


def _yearly_avg(arr):
    return yearly_average(arr, arr[TIME_WEIGHTS_STR])


def _yearly_stdev(arr):
    return _yearly_avg(arr).std(YEAR_STR)


no_reduc = GriddedReduction(_null_reduc, 'full')
time_avg = GriddedReduction(_time_avg, 'av')
yearly_ts = GriddedReduction(_yearly_avg, 'ts')
yearly_stdev = GriddedReduction(_yearly_stdev, 'std')


reg_no_reduc = RegionReduction(_null_reduc, 'full')
reg_time_avg = RegionReduction(_time_avg, 'av')
reg_yearly_ts = RegionReduction(_yearly_avg, 'ts')
reg_yearly_stdev = RegionReduction(_yearly_stdev, 'std',
                                   coord_reduction=reg_yearly_ts)


_TIME_DEFINED_REDUCTIONS = ('av', 'ts', 'std', 'reg.av', 'reg.ts', 'reg.std')
