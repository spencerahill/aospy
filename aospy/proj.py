"""Functionality for representing individual projects for use in aospy."""
import logging
import time


class Proj(object):
    """An object that describes a single project that will use aospy.

    This is the top-level class in the aospy hierarchy of data representations.
    It is meant to contain all of the `Model`, `Run`, and `Region` objects that
    are of relevance to a particular research project.  (Any of these may be
    used by multiple `Proj` objects.)

    The `Proj` class itself provides little functionality, but it is an
    important means of organizing a user's work across different
    projects.  In particular, the output of all calculations using
    `aospy.Calc` are saved in a directory structure whose root is that
    of the `Proj` object specified for that calculation.


    Attributes
    ----------
    name : str
        The run's name
    description : str
        A description of the run
    direc_out, tar_direc_out : str
        The paths to the root directories of, respectively, the standard and
        .tar versions of the output of aospy calculations saved to disk.
    models : dict
        A dictionary with entries of the form ``{model_obj.name: model_obj}``,
        for each of this ``Proj``'s child model objects
    default_models : dict
        The default model objects on which to perform calculations via
        `aospy.Calc` if not otherwise specified
    regions : dict
        A dictionary with entries of the form ``{regin_obj.name: region_obj}``,
        for each of this ``Proj``'s child region objects
"""

    def __init__(self, name, description=None, models=None,
                 default_models=None, regions=None, direc_out='',
                 tar_direc_out=''):
        """
        Parameters
        ----------
        name : str
            The project's name.  This should be unique from that of any other
            `Proj` objects being used.
        description : str, optional
            A description of the model.  This is not used internally by
            aospy; it is solely for the user's information.
        regions : {None, sequence of aospy.Region objects}, optional
            The desired regions over which to perform region-average
            calculations.
        models : {None, sequence of aospy.Model objects}, optional
            The child Model objects of this project.
        default_models : {None, sequence of aospy.Run objects}, optional
            The subset of this Model's runs over which to perform calculations
            by default.
        direc_out, tar_direc_out : str
            Path to the root directories of where, respectively, regular output
            and a .tar-version of the output will be saved to disk.

        Notes
        -----
        Instantiating a `Proj` object has the side-effect of setting the `proj`
        attribute of each of it's child `Model` objects to itself.

        See Also
        --------
        aospy.Model, aospy.Region, aospy.Run

        """
        logging.debug("Initializing Project instance: %s (%s)"
                      % (name, time.ctime()))
        self.name = name
        self.description = '' if description is None else description
        self.direc_out = direc_out
        self.tar_direc_out = tar_direc_out

        if models is None:
            self.models = []
        else:
            self.models = models
        if default_models is None:
            self.default_models = []
        else:
            self.default_models = default_models
        if regions is None:
            self.regions = []
        else:
            self.regions = regions

        # Set the `proj` attribute of the children models.
        for model in self.models:
            setattr(model, 'proj', self)

    def __str__(self):
        return 'Project instance "' + self.name + '"'

    __repr__ = __str__
