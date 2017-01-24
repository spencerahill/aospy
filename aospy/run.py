"""Functionality for storing attributes of a model run or obs product."""
from .utils.times import datetime_or_default


class Run(object):
    """An object that describes a single model 'run' (i.e. simulation).

    Each `Run` object is associated with a parent `Model` object.  This
    parent attribute is not set by `Run` itself, however; it is set during
    the instantation of the parent `Model` object.

    If aospy is being used to work with non climate-model data, the `Run`
    object can be used e.g. to represent different versions of a gridded
    observational data product, with the parent `Model` representing that
    data product more generally.

    Attributes
    ----------
    name : str
        The run's name
    description : str
        A description of the run
    proj : {None, aospy.Proj}
        The run's parent aospy.Proj object
    default_start_date, default_end_date : datetime.datetime
        The default start and end dates of any calculations using this Run
    data_loader : aospy.DataLoader
        The aospy.DataLoader object used to find data on disk corresponding
        to this Run object

    """

    def __init__(self, name=None, description=None, proj=None,
                 default_start_date=None,
                 default_end_date=None,
                 data_loader=None):
        """Instantiate a `Run` object.

        Parameters
        ----------
        name : str
            The run's name.  This must be unique from that of any other
            `Run` objects being used by the parent `Model`.
        description : str, optional
            A description of the model.  This is not used internally by
            aospy; it is solely for the user's information.
        proj : {None, aospy.Proj}, optional
            The parent Proj object.
        data_loader : aospy.DataLoader
            The `DataLoader` object used to find the data on disk to be used
            as inputs for aospy calculations for this run.
        default_start_date, default_end_date : datetime.datetime, optional
            Default start and end dates of calculations to be performed using
            this Model.

        See Also
        --------
        aospy.DataLoader, aospy.Model

        """
        self.name = '' if name is None else name
        self.description = '' if description is None else description
        self.proj = proj

        self.default_start_date = datetime_or_default(
            default_start_date, getattr(data_loader, 'data_start_date', None))
        self.default_end_date = datetime_or_default(
            default_end_date, getattr(data_loader, 'data_end_date', None))

        self.data_loader = data_loader

    def __str__(self):
        return 'aospy.Run instance "%s"' % self.name

    __repr__ = __str__
