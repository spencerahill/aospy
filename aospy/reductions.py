"""Functionality relating to spatiotemporal reduction operations."""
from .internal_names import YEAR_STR


_TIME_DEFINED_REDUCTIONS = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']


registry = {}


def register_reduction(reduction):
    label = reduction.label
    if label in registry:
        raise KeyError("A Reduction object already exists with the desired "
                       "label '{0}', and duplicates are not allowed.  "
                       "The existing object is: "
                       "{1}".format(label, registry[label]))
    else:
        registry[label] = reduction


class RegisterReduction(type):
    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        register_reduction(cls)
        return cls


class Reduction(object, metaclass=RegisterReduction):
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
