import typing as tp

import astropy.io.fits as fits
import astropy.time as astro_time
import astropy.units as u
import h5py
import numpy as np
import scipy.interpolate as interpolate
from astropy.coordinates import get_body, get_moon, get_sun

from gbmgeometry.utils.gbm_time import GBMTime


class PositionInterpolator(object):
    def __init__(self,
                 quats: tp.Iterable[float],
                 sc_pos: tp.Iterable[float],
                 time: float,
                 trigtime=None,
                 factor=1,
                 flags=None):
        """
        An object for interpolating the orbit of 
        GBM

        :param quaternions
        :param sc_pos: 
        :param time:
        :param trigtime:
        :param factor:
        :param flags:
        :returns: 
        :rtype: 

        """

        self._quats: tp.Iterable[float] = quats
        self._sc_pos: tp.Iterable[float] = sc_pos
        self._time: float = time
        self._trigtime = trigtime
        self._factor = factor
        self._flags = flags

        # Interpolate the stuf
        self._interpolate_quaternion()
        self._interpolate_sc_pos()
        self._interpolate_flags()

    @classmethod
    def from_trigdat_hdf5(cls, trigdat_file: str):
        """
        create and interpolator from a trigdat
        HDF5 file

        :param cls: 
        :param trigdat_file: 
        :returns: 
        :rtype: 

        """

        with h5py.File(trigdat_file, "r") as f:

            quats = f["quats"][()]
            sc_pos = f["sc_pos"][()]
            time = f["time"][()]

            trigtime = f.attrs["trigtime"]

        factor = 1

        return cls(
            quats=quats, sc_pos=sc_pos, time=time, trigtime=trigtime, factor=factor
        )

    @classmethod
    def from_poshist_hdf5(cls, poshist_file, T0=None):
        """
        create and interpolator from a poshist
        HDF5 file

        :param cls: 
        :param poshist_file: 
        :param T0:
        :returns: 
        :rtype: 

        """

        with h5py.File(poshist_file, "r") as f:

            quats = f["quats"][()]
            sc_pos = f["sc_pos"][()]
            time = f["time"][()]
            flags = f["flags"][()]

        factor = (u.m).to(u.km)

        if T0 is not None:
            time -= T0

            trigtime = T0

        else:

            trigtime = None

        return cls(
            quats=quats,
            sc_pos=sc_pos,
            time=time,
            trigtime=trigtime,
            factor=factor,
            flags=flags,
        )

    @classmethod
    def from_trigdat(cls, trigdat_file):
        """

        :param cls: 
        :param trigdat_file: 
        :returns: 
        :rtype: 

        """

        with fits.open(trigdat_file) as trigdat:

            trigtime = trigdat["EVNTRATE"].header["TRIGTIME"]
            tstart = trigdat["EVNTRATE"].data["TIME"] - trigtime

            quats = trigdat["EVNTRATE"].data["SCATTITD"]
            sc_pos = trigdat["EVNTRATE"].data["EIC"]

            sort_mask = np.argsort(tstart)
            tstart = tstart[sort_mask]

            quats = quats[sort_mask]
            sc_pos = sc_pos[sort_mask]

            time = tstart

        # the sc is in km so no need to convert

        factor = 1

        return cls(
            quats=quats, sc_pos=sc_pos, time=time, trigtime=trigtime, factor=factor
        )

    @classmethod
    def from_poshist(cls, poshist_file, T0=None):
        """
        create an interpolator from a posthist fits file

        :param cls: 
        :param poshist_file: 
        :param T0: 
        :returns: 
        :rtype: 

        """

        with fits.open(poshist_file) as poshist:

            time = poshist["GLAST POS HIST"].data["SCLK_UTC"]

            flags = poshist["GLAST POS HIST"].data["FLAGS"]

            quats = np.array(
                [
                    poshist["GLAST POS HIST"].data["QSJ_1"],
                    poshist["GLAST POS HIST"].data["QSJ_2"],
                    poshist["GLAST POS HIST"].data["QSJ_3"],
                    poshist["GLAST POS HIST"].data["QSJ_4"],
                ]
            ).T

            sc_pos = np.array(
                [
                    poshist["GLAST POS HIST"].data["POS_X"],
                    poshist["GLAST POS HIST"].data["POS_Y"],
                    poshist["GLAST POS HIST"].data["POS_Z"],
                ]
            ).T

            # if using posthist then units are in m

        factor = (u.m).to(u.km)

        if T0 is not None:
            time -= T0

            trigtime = T0

        else:

            trigtime = None

        return cls(
            quats=quats,
            sc_pos=sc_pos,
            time=time,
            trigtime=trigtime,
            factor=factor,
            flags=flags,
        )

    @property
    def time(self):
        return self._time

    def utc(self, t):

        met = self.met(t)

        time = GBMTime.from_MET(met)
        # print(time.time.fits)
        return time.time.fits

    def astro_time(self, t):

        return astro_time.Time(self.utc(t))

    def met(self, t):

        if self._trigtime is not None:

            met = self._trigtime + t

        else:

            met = t

        return met

    def sun_position(self, t):
        """
        Get the position of the sun at the give time
        relative to the geocenter


        :param t: 
        :returns: 
        :rtype: 

        """

        return get_sun(self.astro_time(t))

    def moon_position(self, t):
        """
        Get the position of the moon at the give time
        relative to the geocenter


        :param t: 
        :returns: 
        :rtype: 

        """

        return get_moon(self.astro_time(t))

    def body_position(self, t, body="uranus"):
        """
        Get the position of the body at the give time
        relative to the geocenter


        :param t: 
        :param body:
        :returns: 
        :rtype: 

        """

        return get_body(body, self.astro_time(t))

    def minmax_time(self):
        """
        return a tuple of (tmin, tmax)
        for the whole interpolator

        :returns: 
        :rtype: 

        """

        return min(self._time), max(self._time)

    def maxtime(self):

        return self._time

    def quaternion(self, t):
        """
        return the quaternion

        :param t: 
        :returns: 
        :rtype: 

        """

        return self._quaternion_t(t).T

    def quaternion_dict(self, t):

        names = [f"quaternion_{i+1}" for i in range(4)]

        return dict(zip(names, self._quaternion_t(t)))

    def sc_pos(self, t):
        """
        the space craft postions in X,Y,Z geocentric coordinates

        units are always assumed to be kilometers

        :param t: 
        :returns: 
        :rtype: 

        """
        return self._scxyz_t(t).T * self._factor

    def sc_pos_dict(self, t):

        names = [f"sc_pos_{s}" for s in ["X", "Y", "Z"]]
        return dict(zip(names, self._scxyz_t(t) * self._factor))

    def is_earth_occulted(self, ra, dec, t):

        earth_radius = 6371.0

        sc_pos = self.sc_pos(t)

        if np.atleast_1d(t).shape[0] == 1:

            fermi_radius = np.sqrt((sc_pos ** 2).sum())

            horizon_angle = 90 - \
                np.rad2deg(np.arccos(earth_radius / fermi_radius))

            min_vis = np.deg2rad(horizon_angle)

            cart_position = ang2cart(ra, dec)

            ang_sep = get_ang(cart_position, -sc_pos)

            return ang_sep < min_vis

        else:

            out = np.ones(len(t), dtype=bool)

            for i, scp in enumerate(sc_pos):

                fermi_radius = np.sqrt((scp ** 2).sum())

                horizon_angle = 90 - \
                    np.rad2deg(np.arccos(earth_radius / fermi_radius))

                min_vis = np.deg2rad(horizon_angle)

                cart_position = ang2cart(ra, dec)

                ang_sep = get_ang(cart_position, -scp)

                out[i] = ang_sep < min_vis

            return out

    def _interpolate_quaternion(self):

        self._quaternion_t = interpolate.interp1d(self._time, self._quats.T)

    def _interpolate_sc_pos(self):

        self._scxyz_t = interpolate.interp1d(self._time, self._sc_pos.T)

    def _interpolate_flags(self):

        if self._flags is not None:

            idx = self._flags == 1

            slices_true = slice_disjoint(idx)

            self._on_times = []

            for s1, s2 in slices_true:

                self._on_times.append([self._time[s1], self._time[s2]])

    def is_fermi_active(self, t):
        """
        test if fermi is active at the time

        :param t: 
        :returns: 
        :rtype: 

        """

        if np.atleast_1d(t).shape[0] == 1:

            if self._flags is None:

                return True

            for tmin, tmax in self._on_times:

                if (tmin <= t) and (t <= tmax):
                    return True

            return False

        else:

            if self._flags is None:

                return np.ones_like(t, dtype=bool)

            out = np.zeros_like(t, dtype=bool)

            for i, tt in enumerate(t):
                for tmin, tmax in self._on_times:

                    if (tmin <= tt) and (tt <= tmax):
                        out[i] = True

            return out

    def sc_matrix(self, t):

        q1, q2, q3, q4 = self.quaternion(t)
        sc_matrix = np.zeros((3, 3))

        sc_matrix[0, 0] = q1 ** 2 - q2 ** 2 - q3 ** 2 + q4 ** 2
        sc_matrix[0, 1] = 2.0 * (q1 * q2 + q4 * q3)
        sc_matrix[0, 2] = 2.0 * (q1 * q3 - q4 * q2)
        sc_matrix[1, 0] = 2.0 * (q1 * q2 - q4 * q3)
        sc_matrix[1, 1] = -(q1 ** 2) + q2 ** 2 - q3 ** 2 + q4 ** 2
        sc_matrix[1, 2] = 2.0 * (q2 * q3 + q4 * q1)
        sc_matrix[2, 0] = 2.0 * (q1 * q3 + q4 * q2)
        sc_matrix[2, 1] = 2.0 * (q2 * q3 - q4 * q1)
        sc_matrix[2, 2] = -(q1 ** 2) - q2 ** 2 + q3 ** 2 + q4 ** 2

        return sc_matrix

    def geo_matrix(self, t):

        return self.sc_matrix(t).T

    def altitude(self, t):
        """

        :param t: 
        :return: 
        """

        earth_radius = 6371.0
        fermi_radius = np.sqrt((self.sc_pos(t) ** 2).sum())

        return fermi_radius - earth_radius

    @staticmethod
    def normalize(x):
        norm = np.sqrt(np.sum(x ** 2, axis=0))
        return x / norm


def ang2cart(ra, dec):
    """
    :param ra:
    :param dec:
    :return:
    """
    pos = np.zeros(3)
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)

    pos[0] = np.cos(dec) * np.cos(ra)
    pos[1] = np.cos(dec) * np.sin(ra)
    pos[2] = np.sin(dec)

    return pos


def get_ang(X1, X2):
    """
    :param X1:
    :param X2:
    :return:
    """
    norm1 = np.sqrt(X1.dot(X1))
    norm2 = np.sqrt(X2.dot(X2))
    tmp = np.clip(np.dot(X1 / norm1, X2 / norm2), -1, 1)

    return np.arccos(tmp)


def slice_disjoint(arr):
    """
    Returns an array of disjoint indices from a bool array
    :param arr: and array of bools
    """

    arr = (arr).nonzero()[0]

    slices = []
    start_slice = arr[0]
    counter = 0
    for i in range(len(arr) - 1):
        if arr[i + 1] > arr[i] + 1:
            end_slice = arr[i]
            slices.append([start_slice, end_slice])
            start_slice = arr[i + 1]
            counter += 1
    if counter == 0:
        return [[arr[0], arr[-1]]]
    if end_slice != arr[-1]:
        slices.append([start_slice, arr[-1]])
    return slices
