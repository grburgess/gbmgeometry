import astropy.time as time
import astropy.units as u
import numpy as np


class GBMTime(object):
    def __init__(self, time_object):

        self._time_object = time_object

        self._current_mjd = self._time_object.mjd

        # get the Fermi MET from the MET

        self._calculate_met_from_mjd()

        self._utc_zero = self._calculate_MJD_from_MET(0)

        # this is when week 9 of the mission starts
        self._utc_start_of_sc_data = "2008-08-07T03:35:44.0"
        self._time_of_start_of_sc_data = time.Time(self._utc_start_of_sc_data)

    @property
    def met(self):

        return self._met

    @property
    def utc(self):

        return self._time_object.iso

    @property
    def time(self):

        return self._time_object

    @property
    def t_zero(self):
        return self._utc_zero

    @property
    def mission_week(self):

        dt = (self._time_object - self._time_of_start_of_sc_data).to(u.week)
        return dt + 10 * u.week

    @classmethod
    def from_UTC_fits(cls, date_string):
        """

        Create a time object from a fits UTC representation


        :param date_string: 
        :return: 
        """

        time_object = time.Time(date_string, format="fits", scale="utc")

        return cls(time_object)

    @classmethod
    def from_MET(cls, met):

        time_object = GBMTime._calculate_MJD_from_MET(met)

        return cls(time_object)

    @staticmethod
    def _calculate_MJD_from_MET(met):

        if met <= 252460801.000:
            utc_tt_diff = 65.184
        elif met <= 362793602.000:
            utc_tt_diff = 66.184
        elif met <= 457401603.000:
            utc_tt_diff = 67.184
        elif met <= 504921604.000:
            utc_tt_diff = 68.184
        else:
            utc_tt_diff = 69.184

        mjdutc = (
            ((met - utc_tt_diff) / 86400.0) + 51910 + 0.0007428703703
        )  # -68.184 added to account for diff between TT and UTC and the 4 leapseconds since 2001
        # mjdtt = ((met) / 86400.0) + 51910 + 0.00074287037037

        return time.Time(mjdutc, scale="utc", format="mjd")

    def _calculate_met_from_mjd(self):
        """
        calculated the Fermi MET given MJD
        :return: 
        """

        if self._current_mjd <= 54832.00000000:
            utc_tt_diff = 65.184
        elif self._current_mjd <= 56109.00000000:
            utc_tt_diff = 66.184
        elif self._current_mjd <= 57204.00000000:
            utc_tt_diff = 67.184
        elif self._current_mjd <= 57754.00000000:
            utc_tt_diff = 68.184
        else:
            utc_tt_diff = 69.184

        self._met = (
            self._current_mjd - 51910 - 0.0007428703703
        ) * 86400.0 + utc_tt_diff  # convert it into MET

    def __add__(self, other):

        if isinstance(other, time.TimeDelta):

            new_time = self._time_object + other

        else:

            # assuming second addition

            dt = time.TimeDelta(other, format="sec")

            new_time = self._time_object + dt

        return GBMTime(new_time)

    def __sub__(self, other):

        if isinstance(other, time.TimeDelta):

            new_time = self._time_object - other

        elif isinstance(other, GBMTime):

            dt = self._time_object - other.time

            return dt

        else:

            # assuming second addition

            dt = time.TimeDelta(other, format="sec")

            new_time = self._time_object - dt

        return GBMTime(new_time)


# def mission_week(met):
