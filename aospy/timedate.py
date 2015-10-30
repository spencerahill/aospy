"""Module for handling times, dates, etc."""
import datetime

import numpy as np
import pandas as pd
import xray


class TimeManager(object):
    """Convert input time specifications into arrays of datetime objects."""

    NUMPYMINYEAR = 1678
    # We start at year 1 (not 0) so this needs to be 1900-1.
    YEAROFFSET = 1899

    @classmethod
    def to_datetime(cls, obj):
        """Create datetime object from the inputted object."""
        if isinstance(obj, bool):
            return obj
        if isinstance(obj, datetime.datetime):
            return obj
        if isinstance(obj, int):
            return datetime.datetime(obj, 1, 1)
        if isinstance(obj, str):
            return cls.str_to_datetime(obj)
        try:
            return datetime.datetime(obj)
        except:
            return cls.str_to_datetime(obj.__str__)

    @staticmethod
    def month_indices(months):
        """Convert string labels for months to integer indices.

        :param months: String matching either 'ann' or some subset of
                       'jfmamjjasond'.  If 'ann', use all months.  Otherwise,
                       use the specified months.
        :type months: str or int
        """
        assert isinstance(months, (int, str))
        if isinstance(months, int):
            return [months]
        if months.lower() == 'ann':
            return range(1, 13)
        first_letter = 'jfmamjjasond'*2
        # Python indexing starts at 0; month indices start at 1 for January.
        st_ind = first_letter.find(months.lower()) + 1
        return range(st_ind, st_ind + len(months))

    def __init__(self, start_date, end_date, months):
        """Instantiate a TimeManager object."""
        self.start_date = self.to_datetime(start_date)
        self.end_date = self.to_datetime(end_date)
        self.months = self.month_indices(months)

    @staticmethod
    def _construct_month_conditional(time, months):
        """Create a conditional statement for selecting data in a DataArray."""
        cond = False
        for month in months:
            cond |= (time['time.month'] == month)
        return cond

    def create_time_array(self):
        """Create an xray.DataArray comprising the desired months."""
        all_months = pd.date_range(start=self.apply_year_offset(self.start_date),
                                   end=self.apply_year_offset(self.end_date), freq='M')
        time = xray.DataArray(all_months, dims=['time'])
        month_cond = self._construct_month_conditional(time, self.months)
        return time[month_cond]

    @staticmethod
    def str_to_datetime(string):
        """Convert a YYYY-MM-DD[...] string to a datetime.datetime object.

        Retains only the year, month, and day from the original object.
        """
        # Retain only the first 10 chars: YYYY-MM-DD
        return datetime.datetime(*[int(d) for d in string[:10].split('-')])

    @staticmethod
    def ymd_to_numpy(year, month, day):
        """Create a numpy.datetime64 with the given year, month, and day."""
        return np.datetime64(
            '{:04d}-{:02d}-{:02d}'.format(year, month, day)
        )

    @classmethod
    def apply_year_offset(cls, date):
        """"Offset the years of a date by default value.

        Hack to address np.datetime64 not supporting years before 1678.
        """
        if date.year < cls.NUMPYMINYEAR:
            offset = cls.YEAROFFSET
        else:
            offset = 0
        return pd.to_datetime(cls.ymd_to_numpy(date.year + offset,
                                               date.month, date.day))


def _get_time(time, start_date, end_date, months, indices=False):
    """Determine indices/values of a time array within the specified interval.

    Assumes time is an xray DataArray and that it can be represented
    by numpy datetime64 objects (i.e. the year is between 1678 and 2262).
    """
    inds = TimeManager._construct_month_conditional(time, months)
    inds &= (time['time'] <= np.datetime64(end_date))
    inds &= (time['time'] >= np.datetime64(start_date))
    if indices == 'only':
        return inds
    elif indices:
        return (inds, time.sel(time=inds))
    else:
        return time.sel(time=inds)
