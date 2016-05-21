"""`Run` class; for storing attributes of a model run or obs product."""

from .timedate import TimeManager


class Run(object):
    """Model run parameters."""
    def _set_direc(self, data_in_direc, ens_mem_prefix, ens_mem_ext,
                   ens_mem_suffix):
        """Set the list of paths containing the Run's netCDF data."""
        if all((ens_mem_prefix, ens_mem_ext, ens_mem_suffix)):
            return [ens_mem_prefix + ext + ens_mem_suffix
                    for ext in ens_mem_ext]
        return data_in_direc

    def __init__(self, name=None, description=None, proj=False,
                 data_in_direc=False, data_in_dur=False,
                 data_in_start_date=False, data_in_end_date=False,
                 data_in_dir_struc='gfdl', data_in_suffix=False,
                 data_in_files=None, default_date_range=False,
                 ens_mem_prefix=False, ens_mem_ext=False, ens_mem_suffix=False,
                 tags=None, idealized=False):
        """Instantiate a `Run` object."""
        self.name = '' if name is None else name
        self.description = '' if description is None else description
        self.proj = proj

        self.data_in_dur = data_in_dur
        self.data_in_dir_struc = data_in_dir_struc
        self.data_in_suffix = data_in_suffix
        self.data_in_files = {} if data_in_files is None else data_in_files
        self.data_in_start_date = TimeManager.to_datetime(data_in_start_date)
        self.data_in_end_date = TimeManager.to_datetime(data_in_end_date,
                                                        end_of_year=True)
        try:
            self.default_date_range = tuple([
                TimeManager.to_datetime(d, end_of_year=eoy)
                for d, eoy in zip(default_date_range, [False, True])]
            )
        except TypeError:
            self.default_date_range = (self.data_in_start_date,
                                       self.data_in_end_date)
        self.tags = [] if tags is None else tags
        self.idealized = idealized
        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        self.data_in_direc = self._set_direc(data_in_direc, ens_mem_prefix,
                                             ens_mem_ext, ens_mem_suffix)

    def __str__(self):
        return 'Run instance "%s"' % self.name

    __repr__ = __str__
