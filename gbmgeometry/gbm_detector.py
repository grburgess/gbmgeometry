__author__ = "drjfunk"
import numpy as np
from astropy.coordinates import (
    SkyCoord,
    get_sun,
    get_body,
    get_body_barycentric,
    get_moon,
)
from spherical_geometry.polygon import SphericalPolygon
import astropy.units as u

from gbmgeometry.utils.plotting import skyplot, SphericalCircle

from .gbm_frame import GBMFrame


class GBMDetector(object):
    def __init__(self, name, quaternion, sc_pos=None, time=None):
        """
        
        :param name: 
        :param quaternion: 
        :param sc_pos: 
        :param time: 
        """

        self._name = name

        self.update_position(quaternion, sc_pos, time)

        self._xyz = np.array(
            [
                np.cos(self._az) * np.sin(self._zen),
                np.sin(self._az) * np.sin(self._zen),
                np.cos(self._zen),
            ]
        )

    def update_position(self, quaternion, sc_pos=None, time=None):
        """
        
        :param quaternion: 
        :param sc_pos: 
        :param time: 
        :return: 
        """

        self._time = time

        q1, q2, q3, q4 = quaternion

        if sc_pos is not None:
            scx, scy, scz = sc_pos
            self._heigth = np.sqrt(scx * scx + scy * scy + scz * scz)

        else:
            scx = None
            scy = None
            scz = None
        self._quaternion = quaternion
        self._sc_pos = sc_pos
        # define the direction of the detector in the GBMFrame
        self._center = SkyCoord(
            lon=self._az * u.deg,
            lat=self._zen * u.deg,
            unit="deg",
            frame=GBMFrame(
                quaternion_1=q1,
                quaternion_2=q2,
                quaternion_3=q3,
                quaternion_4=q4,
                sc_pos_X=scx,
                sc_pos_Y=scy,
                sc_pos_Z=scz,
            ),
        )

        self._center_icrs = self._center.icrs

        if self._time is not None:
            # we can calculate the sun position
            # in GCRS
            tmp_sun = get_sun(self._time)
            # in ICRS
            self._sun_position_icrs = SkyCoord(
                tmp_sun.ra.deg,
                tmp_sun.dec.deg,
                unit="deg",
                frame="gcrs",
                obstime=self._time,
            ).icrs
            # in CenterFrame
            self._sun_position = self._sun_position_icrs.transform_to(
                self._center.frame
            )

            # position of earth in satellite frame:

        if sc_pos is not None:

            self._earth_pos_norm = self.geo_to_gbm(
                -self._sc_pos / np.linalg.norm(self._sc_pos)
            )
            scxn, scyn, sczn = self._earth_pos_norm
            earth_theta = np.arccos(
                sczn / np.sqrt(scxn * scxn + scyn * scyn + sczn * sczn)
            )
            earth_phi = np.arctan2(scyn, scxn)
            earth_ra = np.rad2deg(earth_phi)
            if earth_ra < 0:
                earth_ra = earth_ra + 360
            earth_dec = 90 - np.rad2deg(earth_theta)

            # earth as SkyCoord
            self._earth_position = SkyCoord(
                lon=earth_ra * u.deg,
                lat=earth_dec * u.deg,
                unit="deg",
                frame=GBMFrame(
                    quaternion_1=q1,
                    quaternion_2=q2,
                    quaternion_3=q3,
                    quaternion_4=q4,
                    sc_pos_X=scx,
                    sc_pos_Y=scy,
                    sc_pos_Z=scz,
                ),
            )
        self._scx = scx
        self._scy = scy
        self._scz = scz

    def set_quaternion(self, quaternion):
        """
        Parameters
        ----------
        quaternion

        """

        self._quaternion = quaternion

        q1, q2, q3, q4 = quaternion

        if self._sc_pos is not None:
            scx, scy, scz = self._sc_pos

        else:
            scx = None
            scy = None
            scz = None

        self._center = SkyCoord(
            lon=self._az * u.deg,
            lat=self._zen * u.deg,
            unit="deg",
            frame=GBMFrame(
                quaternion_1=q1,
                quaternion_2=q2,
                quaternion_3=q3,
                quaternion_4=q4,
                sc_pos_X=scx,
                sc_pos_Y=scy,
                sc_pos_Z=scz,
            ),
        )

        self._center_icrs = self._center.icrs

    def set_sc_pos(self, sc_pos):
        """
        Parameters
        ----------
        quaternion

        """

        q1, q2, q3, q4 = self._quaternion

        if sc_pos is not None:
            scx, scy, scz = sc_pos

        else:
            scx = None
            scy = None
            scz = None

        self._center = SkyCoord(
            lon=self._az * u.deg,
            lat=self._zen * u.deg,
            unit="deg",
            frame=GBMFrame(
                quaternion_1=q1,
                quaternion_2=q2,
                quaternion_3=q3,
                quaternion_4=q4,
                sc_pos_X=scx,
                sc_pos_Y=scy,
                sc_pos_Z=scz,
            ),
        )

        self._center_icrs = self._center.icrs

    def get_fov(self, radius, fermi_frame=False):
        """
        Returns
        -------
        array of RA and DEC

        """

        steps = 500

        if fermi_frame:
            fermi = self._center

            poly = SphericalPolygon.from_cone(
                fermi.lon.value, fermi.lat.value, radius, steps=steps
            )

        else:

            j2000 = self._center.icrs

            poly = SphericalPolygon.from_cone(
                j2000.ra.value, j2000.dec.value, radius, steps=steps
            )

        # ra, dec
        return [p for p in poly.to_radec()][0]

    def get_center(self):
        return self._center

    @property
    def height(self):  # heigth above earth center
        """
        :return: height of the satelite
        """
        return self._heigth

    @property
    def sun_position(self):
        """
        :return: sun position as SkyCoord object in sat frame
        """

        return self._sun_position

    @property
    def sun_position_icrs(self):
        """
        :return: sun position as SkyCood object in icrs frame
        """
        return self._sun_position.icrs

    @property
    def moon_position_icrs(self):
        """                                                                                                                                                                                                                                                                   
        :return: moon position as SkyCood object in icrs frame                                                                                                                                                                                                                
        """
        tmp_moon = get_moon(self._time)
        # in ICRS
        return SkyCoord(
            tmp_moon.ra.deg,
            tmp_moon.dec.deg,
            unit="deg",
            frame="gcrs",
            obstime=self._time,
        ).icrs

    @property
    def sun_angle(self):
        """
        :return: angle between det pointing and sun position
        """

        return self._center.separation(self._sun_position)

    @property
    def sun_earth_angle(self):
        """
        :return: angle between center of earth an sun postition (seen from the satelite)
        """
        return self._sun_position.separation(self._earth_position)

    @property
    def earth_position(self):
        """
        :return: earth position as SkyCoord in sat frame
        """

        return self._earth_position

    @property
    def earth_position_icrs(self):
        """
        :return: earth position as SkyCoord in ICRS frame
        """
        return self._earth_position.icrs

    @property
    def earth_az_zen_sat(self):
        """
        :return: the longitude and latitude of the earth center in sat. frame
        """
        return [self._earth_position.lon.deg, self._earth_position.lat.deg]

    @property
    def det_ra_dec_icrs(self):
        """
        :return: the det pointing in icrs frame (ra, dec
        """
        return [(self._center.icrs).ra.deg, (self._center.icrs).dec.deg]

    def geo_to_gbm(self, pos_geo):
        """ Compute the transformation from heliocentric Sgr coordinates to
            spherical Galactic.
        """
        q1, q2, q3, q4 = self._quaternion  # q1,q2,q3,q4
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

        X0 = np.dot(sc_matrix[0, :], pos_geo)
        X1 = np.dot(sc_matrix[1, :], pos_geo)
        X2 = np.clip(np.dot(sc_matrix[2, :], pos_geo), -1.0, 1.0)
        pos_sat = [X0, X1, X2]
        return np.array(pos_sat)

    def gbm_to_geo(self, pos_gbm):
        """ Compute the transformation from heliocentric Sgr coordinates to
            spherical Galactic.
        """
        q1, q2, q3, q4 = self._quaternion  # q1,q2,q3,q4
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

        X0 = np.dot(sc_matrix[:, 0], pos_gbm)
        X1 = np.dot(sc_matrix[:, 1], pos_gbm)
        X2 = np.clip(np.dot(sc_matrix[:, 2], pos_gbm), -1.0, 1.0)
        pos_geo = [X0, X1, X2]
        return np.array(pos_geo)

    @property
    def earth_angle(self):
        """
        :return: angle between det pointing and earth position
        """
        return self._center.separation(self._earth_position)

    @property
    def sun_pos(self):
        """
        :return: sun position in sat frame expressed with x,y and z coordinate
        """
        lon, lat = self._sun_position.lon.deg, self._sun_position.lat.deg
        x = np.cos(lon) * np.cos(lat)
        y = np.sin(lon) * np.cos(lat)
        z = np.sin(lat)
        return np.array([x, y, z])

    @property
    def sun_lon_lat(self):
        """                                                                                                                                                                                                                                                                   
        :return: sun position in sat frame expressed with lon, lat                                                                                                                                                                                                 
        """
        return self._sun_position.lon.deg, self._sun_position.lat.deg

    @property
    def earth_pos(self):
        """                                                                                                                                                                                                                                                                   
        :return: earth position in sat frame expressed with x,y and z coordinate                                                                                                                                                                                             
        """
        lon, lat = self._earth_position.lon.deg, self._earth_position.lat.deg
        x = np.cos(lon) * np.cos(lat)
        y = np.sin(lon) * np.cos(lat)
        z = np.sin(lat)
        return np.array([x, y, z])

    @property
    def center(self):

        return self._center

    @property
    def center_icrs(self):

        return self._center_icrs

    @property
    def name(self):

        return self._name

    @property
    def az(self):

        return self._az

    @property
    def zen(self):

        return self._zen

    @property
    def xyz(self):

        return self._xyz

    @property
    def mount_point(self):

        return self._mount_point

    def plot_pointing(self, ax=None, fov=None, **kwargs):

        if ax is None:

            skymap_kwargs = {}

            if "projection" in kwargs:
                skymap_kwargs["projection"] = kwargs.pop("projection")

            if "center" in kwargs:
                skymap_kwargs["center"] = kwargs.pop("center")

            if "radius" in kwargs:
                skymap_kwargs["radius"] = kwargs.pop("radius")

            ax = skyplot(**skymap_kwargs)

        # compute the annulus for this set of detectors

        if fov is not None:

            assert fov > 0, "fov must be a positive number in units of deg"

            circle = SphericalCircle(
                self._center_icrs.ra,
                self._center_icrs.dec,
                fov,
                vertex_unit=u.deg,
                resolution=5000,
                #            edgecolor=color,
                transform=ax.get_transform("icrs"),
                **kwargs
            )

            ax.add_patch(circle)

        else:

            ax.scatter(
                self._center_icrs.ra.deg,
                self.center_icrs.dec.deg,
                transform=ax.get_transform("icrs"),
                **kwargs
            )

        return ax


class NaI0(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 45.89
        self._zen = 90 - 20.58

        self._mount_point = np.array([96.1, 80.4, 107.6])

        super(NaI0, self).__init__("n0", quaternion, sc_pos, time)


class NaI1(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 45.11
        self._zen = 90 - 45.31
        self._mount_point = np.array([101.1, 72.8, 72.1])

        super(NaI1, self).__init__("n1", quaternion, sc_pos, time)


class NaI2(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 58.44
        self._zen = 90 - 90.21
        self._mount_point = np.array([109.0, 58.1, 99.0])

        super(NaI2, self).__init__("n2", quaternion, sc_pos, time)


class NaI3(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 314.87
        self._zen = 90 - 45.24
        self._mount_point = np.array([97.7, -76.3, 102.5])

        super(NaI3, self).__init__("n3", quaternion, sc_pos, time)


class NaI4(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 303.15
        self._zen = 90.0 - 90.27
        self._mount_point = np.array([109.0, -57.5, 83.6])

        super(NaI4, self).__init__("n4", quaternion, sc_pos, time)


class NaI5(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 3.35
        self._zen = 90 - 89.97
        self._mount_point = np.array([99.6, -49.7, 100.1])

        super(NaI5, self).__init__("n5", quaternion, sc_pos, time)


class NaI6(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 224.93
        self._zen = 90 - 20.43
        self._mount_point = np.array([-95.8, -80.3, 107.1])

        super(NaI6, self).__init__("n6", quaternion, sc_pos, time)


class NaI7(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 224.62
        self._zen = 90 - 46.18
        self._mount_point = np.array([-100.6, -72.5, 71.6])

        super(NaI7, self).__init__("n7", quaternion, sc_pos, time)


class NaI8(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 236.61
        self._zen = 90 - 89.97
        self._mount_point = np.array([-108.4, -57.2, 99.0])

        super(NaI8, self).__init__("n8", quaternion, sc_pos, time)


class NaI9(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 135.19
        self._zen = 90 - 45.55
        self._mount_point = np.array([-97.5, 76.5, 102.5])

        super(NaI9, self).__init__("n9", quaternion, sc_pos, time)


class NaIA(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """

        self._az = 123.73
        self._zen = 90 - 90.42
        self._mount_point = np.array([-108.7, 57.7, 83.7])

        super(NaIA, self).__init__("na", quaternion, sc_pos, time)


class NaIB(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 183.74
        self._zen = 90 - 90.32
        self._mount_point = np.array([-99.3, 50.0, 100.2])

        super(NaIB, self).__init__("nb", quaternion, sc_pos, time)


class BGO0(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 0.0
        self._zen = 0.0
        self._mount_point = np.array([126.05, 0.13, 63.32])

        super(BGO0, self).__init__("b0", quaternion, sc_pos, time)


class BGO1(GBMDetector):
    def __init__(self, quaternion, sc_pos=None, time=None):
        """

        Parameters
        ----------
        quaternion
        """
        self._az = 180.0
        self._zen = 0.0
        self._mount_point = np.array([-126.14, 0.01, 67.22])

        super(BGO1, self).__init__("b1", quaternion, sc_pos, time)
