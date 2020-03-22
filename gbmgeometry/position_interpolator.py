import astropy.io.fits as fits
import h5py
import astropy.units as u
import numpy as np
import scipy.interpolate as interpolate

from gbmgeometry.utils.gbm_time import GBMTime


class PositionInterpolator(object):
    def __init__(self, quats, sc_pos, time, trigtime=None, factor=1):
        """FIXME! briefly describe function

        :param poshist: 
        :param T0: 
        :param trigdat: 
        :returns: 
        :rtype: 

        """

        self._quats = quats
        self._sc_pos = sc_pos
        self._time = time
        self._trigtime = trigtime
        self._factor = factor

        # Interpolate the stuf
        self._interpolate_quaternion()
        self._interpolate_sc_pos()

    @classmethod
    def from_trigdat_hdf5(cls, trigdat_file):
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

        factor = (u.m).to(u.km)

        if T0 is not None:
            time -= T0

            trigtime = T0

        else:

            trigtime = None

        return cls(
            quats=quats, sc_pos=sc_pos, time=time, trigtime=trigtime, factor=factor
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
            quats=quats, sc_pos=sc_pos, time=time, trigtime=trigtime, factor=factor
        )

    def utc(self, t):

        if self._trigtime is not None:

            met = self._trigtime + t

        else:

            met = t

        time = GBMTime.from_MET(met)
        # print(time.time.fits)
        return time.time.fits

    def met(self, t):

        if self._trigtime is not None:

            met = self._trigtime + t

        else:

            met = t

        return met

    def maxtime(self):

        return self._time

    def quaternion(self, t):
        """
        Gets an interpolated quaternion as a function of time

        Parameters
        ----------
        t

        Returns
        -------
        A Fermi GBM quaternion

        """

        return self._quaternion_t(t)

    def quaternion_dict(self, t):

        names = [f"quaternion_{i+1}" for i in range(4)]

        return dict(zip(names, self._quaternion_t(t)))

    def sc_pos(self, t):
        """

        Parameters
        ----------
        t

        Returns
        -------
        Fermi GBM spacecraft position

        """
        return self._scxyz_t(t) * self._factor

    def sc_pos_dict(self, t):

        names = [f"sc_pos_{s}" for s in ["X", "Y", "Z"]]
        return dict(zip(names, self._scxyz_t(t) * self._factor))

    def _interpolate_quaternion(self):

        self._quaternion_t = interpolate.interp1d(self._time, self._quats.T)

    def _interpolate_sc_pos(self):

        self._scxyz_t = interpolate.interp1d(self._time, self._sc_pos.T)

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
