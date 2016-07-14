__author__ = "drjfunk"
from .gbm_frame import GBMFrame
from astropy.coordinates import SkyCoord
from spherical_geometry.polygon import SphericalPolygon


class GBMDetector(object):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._center = SkyCoord(Az=self._az, Zen=self._zen, unit='deg', frame=GBMFrame(quaternion=quaternion))

        self._quaternion = quaternion



    def set_quaternion(self, quaternion):
        """
        Parameters
        ----------
        quaternion

        """

        self._quaternion = quaternion

        self._center = SkyCoord(Az=self._az, Zen=self._zen, unit='deg', frame=GBMFrame(quaternion=quaternion))

    def get_fov(self, radius):
        """
        Returns
        -------
        array of RA and DEC

        """

        steps = 500

        j2000 = self._center.icrs

        poly = SphericalPolygon.from_cone(j2000.ra.value,
                                          j2000.dec.value,
                                          radius,
                                          steps=steps)

        # ra, dec
        return [p for p in poly.to_radec()][0]

    def get_center(self):
        return self._center


class NaI0(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 45.89
        self._zen = 20.58 - 90.

        super(NaI0, self).__init__(quaternion)


class NaI1(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 45.11
        self._zen = 45.31 - 90.

        super(NaI1, self).__init__(quaternion)


class NaI2(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 58.44
        self._zen = 90.21 - 90.

        super(NaI2, self).__init__(quaternion)


class NaI3(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 314.87
        self._zen = 45.24 - 90.

        super(NaI3, self).__init__(quaternion)


class NaI4(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 303.15
        self._zen = 90.27 - 90.

        super(NaI4, self).__init__(quaternion)


class NaI5(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 3.35
        self._zen = 89.97 - 90.

        super(NaI5, self).__init__(quaternion)


class NaI6(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 224.93
        self._zen = 20.43 - 90.

        super(NaI6, self).__init__(quaternion)


class NaI7(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 224.62
        self._zen = 46.18 - 90.

        super(NaI7, self).__init__(quaternion)


class NaI8(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 236.61
        self._zen = 89.97 - 90.

        super(NaI8, self).__init__(quaternion)


class NaI9(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 135.19
        self._zen = 45.55 - 90.

        super(NaI9, self).__init__(quaternion)


class NaIA(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """

        self._az = 123.73
        self._zen = 90.42 - 90.

        super(NaIA, self).__init__(quaternion)


class NaIB(GBMDetector):
    def __init__(self, quaternion):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 183.74
        self._zen = 90.32 - 90.

        super(NaIB, self).__init__(quaternion)
