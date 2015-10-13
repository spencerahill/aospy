"""Module for handling times, dates, etc."""
import numpy as np
import pandas as pd
import xray


class TimeManager(object):
    """Convert input time specifications into arrays of datetime objects."""

    datetime_type = np.datetime64

    @staticmethod
    def to_datetime(obj, datetime_type):
        """Create datetime object from the inputted object."""
        try:
            return datetime_type(obj)
        except:
            raise

    def __init__(self, start_date, end_date, months):
        """Instantiate a TimeManager object."""
        self.start_date = self.to_datetime(start_date, self.datetime_type)
        self.end_date = self.to_datetime(end_date, self.datetime_type)
        self.months = months

    def _construct_month_conditional(self, time, months):
        """Create a conditional statement for selecting data in a DataArray."""
        cond = False
        for month in months:
            cond |= (time['time.month'] == month)
        return cond

    def create_time_array(self, start_date, end_date, months):
        """Create an xray.DataArray comprising the desired months."""
        all_months = pd.date_range(start=start_date, end=end_date, freq='M')
        time = xray.DataArray(all_months, dims=['time'])
        month_conditional = self._construct_month_conditional(time, months)
        return time[month_conditional]
