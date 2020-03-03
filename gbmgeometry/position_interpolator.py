import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import scipy.interpolate as interpolate

from gbmgeometry.utils.gbm_time import GBMTime


class PositionInterpolator(object):
    def __init__(self, poshist=None, T0=None, trigdat=None):

        # Use position history file
        """

        Parameters
        ----------
        poshist
        T0
        trigdat

        Returns
        -------

        """
        if poshist is not None:

            with fits.open(poshist) as poshist:
                # poshist = fits.open(poshist)
                self._time = poshist["GLAST POS HIST"].data["SCLK_UTC"]
                self._quats = np.array(
                    [
                        poshist["GLAST POS HIST"].data["QSJ_1"],
                        poshist["GLAST POS HIST"].data["QSJ_2"],
                        poshist["GLAST POS HIST"].data["QSJ_3"],
                        poshist["GLAST POS HIST"].data["QSJ_4"],
                    ]
                ).T

                self._sc_pos = np.array(
                    [
                        poshist["GLAST POS HIST"].data["POS_X"],
                        poshist["GLAST POS HIST"].data["POS_Y"],
                        poshist["GLAST POS HIST"].data["POS_Z"],
                    ]
                ).T

                # if using posthist then units are in m

                self._factor = (u.m).to(u.km)

                if T0 is not None:
                    self._time -= T0

                    self._trigtime = T0

                else:

                    self._trigtime = None

            # poshist.close()

        elif trigdat is not None:

            with fits.open(trigdat) as trigdat:
                # trigdat = fits.open(trigdat)
                trigtime = trigdat["EVNTRATE"].header["TRIGTIME"]
                tstart = trigdat["EVNTRATE"].data["TIME"] - trigtime

                self._trigtime = trigtime

                self._quats = trigdat["EVNTRATE"].data["SCATTITD"]
                self._sc_pos = trigdat["EVNTRATE"].data["EIC"]

                sort_mask = np.argsort(tstart)
                tstart = tstart[sort_mask]

                self._quats = self._quats[sort_mask]
                self._sc_pos = self._sc_pos[sort_mask]

                self._time = tstart

            # trigdat.close()

            # the sc is in km so no need to convert

            self._factor = 1

        else:

            print("No file passed. Exiting")
            return

        # Interpolate the stuf
        self._interpolate_quaternion()
        self._interpolate_sc_pos()

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
