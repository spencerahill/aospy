"""Module for handling times, dates, etc."""
import datetime

import netCDF4
import numpy as np
import pandas as pd
import xray


class TimeManager(object):
    """Convert input time specifications into arrays of datetime objects."""

    NUMPYMINYEAR = 1678
    YEAROFFSET = 1900

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

    def _construct_month_conditional(self, time, months):
        """Create a conditional statement for selecting data in a DataArray."""
        cond = False
        for month in months:
            cond |= (time['time.month'] == month)
        return cond

    def create_time_array(self):
        """Create an xray.DataArray comprising the desired months."""
        all_months = pd.date_range(start=self.start_date,
                                   end=self.end_date, freq='M')
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


def _get_time_xray(time, start_date, end_date, months, indices=False):
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


def _get_time(time, units, calendar, start_yr, end_yr, months, indices=False):
    """Determine the indices of a time array falling in a specified interval.

    Given a start year, end year, and subset of each year, determine which of
    the input time array's values fall within that range.

    :param time: netCDF4 variable object specifying time
    :param start_yr, end_yr: Start and end years, inclusive, of desired time
                             range.
    :type start_yr, end_yr: int
    :param months: Subset of the annual cycle to sample.
    :type months: Iterable of ints in the range (1,13), inclusive.
    :param indices: Return an array of indices if True, otherwise return
                    the time array itself at those time indices.
    :type indices: bool
    """
    dates = netCDF4.num2date(time[:], units, calendar.lower())
    inds = [i for i, date in enumerate(dates) if (date.month in months) and
            (date.year in range(start_yr, end_yr+1))]
    if indices == 'only':
        return inds
    elif indices:
        return inds, time[inds]
    else:
        return time[inds]
