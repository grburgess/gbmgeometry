import numpy as np
import scipy.interpolate as interpolate
import astropy.io.fits as fits


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

            poshist = fits.open(poshist)
            self._time = poshist['GLAST POS HIST'].data['SCLK_UTC']
            self._quats = np.array([poshist['GLAST POS HIST'].data['QSJ_1'],
                                    poshist['GLAST POS HIST'].data['QSJ_2'],
                                    poshist['GLAST POS HIST'].data['QSJ_3'],
                                    poshist['GLAST POS HIST'].data['QSJ_4']])

            self._sc_pos = np.array([poshist['GLAST POS HIST'].data['POS_X'],
                                     poshist['GLAST POS HIST'].data['POS_Y'],
                                     poshist['GLAST POS HIST'].data['POS_Z']])

            if T0 is not None:
                self._time -= T0

            poshist.close()

        elif trigdat is not None:

            trigdat = fits.open(trigdat)
            trigtime = trigdat['EVNTRATE'].header['TRIGTIME']
            tstart = trigdat['EVNTRATE'].data['TIME'] - trigtime

            self._quats = trigdat['EVNTRATE'].data['SCATTITD']
            self._sc_pos = trigdat['EVNTRATE'].data['EIC']

            sort_mask = np.argsort(tstart)
            tstart = tstart[sort_mask]

            self._quats = self._quats[sort_mask]
            self._sc_pos = self._sc_pos[sort_mask]

            self._time = tstart


            trigdat.close()

        else:

            print "No file passed. Exiting"
            return

        # Interpolate the stuf
        self._interpolate_quaternion()
        self._interpolate_sc_pos()

    def quaternion(self, t):
        """
        Gets an itnerpolated quaternion as a function of time

        Parameters
        ----------
        t

        Returns
        -------
        A Fermi GBM quarternion

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
        return self._scxyz_t(t)

    def _interpolate_quaternion(self):

        self._quaternion_t = interpolate.interp1d(self._time, self._quats.T)

    def _interpolate_sc_pos(self):

        self._scxyz_t = interpolate.interp1d(self._time, self._sc_pos.T)

    @staticmethod
    def normalize(x):
        norm = np.sqrt(np.sum(x ** 2, axis=0))
        return x / norm
