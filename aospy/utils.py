"""aospy.utils: utility functions for the aospy module."""
import numpy as np

from . import user_path

def _load_user_data(name):
    """Load user data from aospy_path for given module name.

    File must be located in the `aospy_path` directory and be the same name
    as the desired aospy module subpackage, namely one of `regions`, `calcs`,
    `variables`, and `projects`.
    """
    import imp
    return imp.load_source(
        name, '/'.join([user_path, name, '__init__.py']).replace('//','/')
    )


def _set_attr_from_init_kwargs(obj, kwargs):
    """
    Given a dict of keyword arguments and their values as input for an
    __init__() call, set each attribute of the object with the name and value
    of the inputted arguments.
    """
    for kwarg_key, kwarg_value in kwargs.iteritems():
        setattr(obj, kwarg_key, kwarg_value)


def _get_parent_attr(obj, attr):
    """
    Check if the object has the given attribute and it is non-empty.  If not,
    check each parent object for the attribute and use the first one found.
    """
    # Get the attribute from the object itself, if it exists.
    try:
        return getattr(obj, attr)
    # Otherwise, loop through the parent objects until finding the attr.
    except AttributeError:
        for parent_obj in ('var', 'run', 'model', 'proj'):
            try:
                return getattr(getattr(obj, parent_obj), attr)
            except AttributeError:
                pass


def _set_named_attr_dict(obj, attr, attr_list):
    """Set attribute lists that are dicts of {name: value} type."""
    setattr(obj, attr, {item.name: item for item in attr_list})


def _set_parent_object_of_child(parent_obj, parent_attr, child):
    setattr(child, parent_attr, parent_obj)


def level_thickness(levs):
    """
    Calculates the thickness, in Pa, of each pressure level.

    Assumes that the pressure values given are at the center of that model
    level, except for the lowest value (typically 1000 hPa), which is the
    bottom boundary. The uppermost level extends to 0 hPa.

    """
    # Bottom level extends from levs[0] to halfway betwen levs[0]
    # and levs[1].
    levs = to_pascal(levs)
    lev_thick = [0.5*(levs[0] - levs[1])]
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    for k in range(1, levs.size-1):
        lev_thick.append(0.5*(levs[k-1] - levs[k+1]))
    # Top level extends from halfway between top two levels to 0 hPa.
    lev_thick.append(0.5*(levs[-2] + levs[-1]))
    # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
    return np.array(lev_thick)
