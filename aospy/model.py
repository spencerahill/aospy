class Model(object):
    """Model and parameters of local data stored and desired."""
    def __init__(self, name='', proj=False, nc_grid_paths=(), nc_dur=False,
                 nc_start_yr=False, nc_end_yr=False, default_yr_range=False,
                 runs=(), default_runs=()):
        self.name = name
        self.description = description
        self.proj = proj
        self.nc_dur = nc_dur
        self.nc_grid_paths = nc_grid_paths
        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.nc_dir_struc = nc_dir_struc
        self.default_yr_range = default_yr_range
        self.runs = runs
        self.default_runs = default_runs

        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        
        # Use the inputted names and netCDF filepath to create grid data.
        _set_mult_nc_grid_attr(self)
        self._set_sfc_area()
        # Populate the runs dictionary with the specified list of runs.
        if 'runs' in kwargs.keys():
            _set_named_attr_dict(self, 'runs', kwargs['runs'])
            # Create dict of default runs to use.
            _set_named_attr_dict(
                self, 'default_runs',
                kwargs.get('default_runs', kwargs['runs'])
                )

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    def _set_levs_thick(self):
        """
        Set the 1D array holding the pressure thickness of the model levels.
        """
        if self.level:
            self.levs_thick = calc_levs_thick(self.level)
        else:
            self.levs_thick = None

    def _set_sfc_area(self):
        """Set the 2D array holding the surface area of gridboxes."""
        if getattr(self, 'sfc_area', None) is not None:
            return
        else:
            try:
                sfc_area = calc_grid_sfc_area(
                    self.lon_bounds, self.lat_bounds, gridcenter=False
                )
            except AttributeError:
                sfc_area = calc_grid_sfc_area(
                    self.lon, self.lat, gridcenter=True
                )
            self.sfc_area = sfc_area

    def _get_dt(self):
        """Get the model's time spacing."""
        import numpy as np
        if self.time_dur is not None:
            return self.time_dur
        elif self.time_bounds is not None:
            return self.time_bounds[:,1] - self.time_bounds[:,0]
        else:
            return np.ones(np.size(self.time))

    def _set_runs(self, runs_list):
        """Create the Model's dictionary of Runs."""
        self.runs = {run.name: run for run in runs_list}
        for run in runs_list:
            run.model = self

def model_inst(model, parent_proj=False):
    """Convert string of an aospy.model name to an aospy.model instance."""
    if parent_proj and type(parent_proj) is not Proj:
        parent_proj = _proj_inst(parent_proj)
    if type(model) is Model:
        model_inst =  model
    elif type(model) is str:
        model_inst = parent_proj.models[model]
    elif type(model) in (list, tuple):
        model_inst = [_model_inst(mod, proj) for (mod, proj)
                      in zip(model, parent_proj)]
        if type(model) is tuple:
            model_inst = tuple(model_inst)
    if parent_proj:
        parent_proj = _proj_inst(parent_proj)
        try:
            model_inst.proj = parent_proj
        except AttributeError:
            pass
    return model_inst
