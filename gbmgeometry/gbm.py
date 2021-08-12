from collections import OrderedDict

import astropy.coordinates as coord
import astropy.units as u

# import mpl_toolkits.basemap as bm
import numpy as np
import pandas as pd
import spherical_geometry.polygon as sp

from gbmgeometry.position_interpolator import PositionInterpolator
from gbmgeometry.utils.gbm_time import GBMTime
from gbmgeometry.utils.plotting import SphericalCircle, skyplot

from .gbm_detector import (
    BGO0,
    BGO1,
    NaI0,
    NaI1,
    NaI2,
    NaI3,
    NaI4,
    NaI5,
    NaI6,
    NaI7,
    NaI8,
    NaI9,
    NaIA,
    NaIB,
)
from .gbm_frame import GBMFrame

# import seaborn as sns

_det_color_cycle = np.linspace(0, 1, 12)


class GBM(object):
    def __init__(self, quaternion, sc_pos=None, gbm_time=None):
        """


        :param quaternion: 
        :param sc_pos: 
        :param gbm_time: 
        :returns: 
        :rtype: 

        """

        if gbm_time is not None:

            if isinstance(gbm_time, str):

                self._gbm_time = GBMTime.from_UTC_fits(gbm_time)

            else:

                # assuming MET

                self._gbm_time = GBMTime.from_MET(gbm_time)

        else:

            self._gbm_time = None

        self._quaternion = quaternion
        self._sc_pos = sc_pos

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
            self._earth_position = coord.SkyCoord(
                lon=earth_ra * u.deg,
                lat=earth_dec * u.deg,
                unit="deg",
                frame=GBMFrame(
                    quaternion_1=quaternion[0],
                    quaternion_2=quaternion[1],
                    quaternion_3=quaternion[2],
                    quaternion_4=quaternion[3],
                    sc_pos_X=sc_pos[0],
                    sc_pos_Y=sc_pos[1],
                    sc_pos_Z=sc_pos[2],
                ),
            )

        if self._gbm_time is not None:

            self.n0 = NaI0(quaternion, sc_pos, self._gbm_time.time)
            self.n1 = NaI1(quaternion, sc_pos, self._gbm_time.time)
            self.n2 = NaI2(quaternion, sc_pos, self._gbm_time.time)
            self.n3 = NaI3(quaternion, sc_pos, self._gbm_time.time)
            self.n4 = NaI4(quaternion, sc_pos, self._gbm_time.time)
            self.n5 = NaI5(quaternion, sc_pos, self._gbm_time.time)
            self.n6 = NaI6(quaternion, sc_pos, self._gbm_time.time)
            self.n7 = NaI7(quaternion, sc_pos, self._gbm_time.time)
            self.n8 = NaI8(quaternion, sc_pos, self._gbm_time.time)
            self.n9 = NaI9(quaternion, sc_pos, self._gbm_time.time)
            self.na = NaIA(quaternion, sc_pos, self._gbm_time.time)
            self.nb = NaIB(quaternion, sc_pos, self._gbm_time.time)
            self.b0 = BGO0(quaternion, sc_pos, self._gbm_time.time)
            self.b1 = BGO1(quaternion, sc_pos, self._gbm_time.time)

        else:

            self.n0 = NaI0(quaternion, sc_pos, None)
            self.n1 = NaI1(quaternion, sc_pos, None)
            self.n2 = NaI2(quaternion, sc_pos, None)
            self.n3 = NaI3(quaternion, sc_pos, None)
            self.n4 = NaI4(quaternion, sc_pos, None)
            self.n5 = NaI5(quaternion, sc_pos, None)
            self.n6 = NaI6(quaternion, sc_pos, None)
            self.n7 = NaI7(quaternion, sc_pos, None)
            self.n8 = NaI8(quaternion, sc_pos, None)
            self.n9 = NaI9(quaternion, sc_pos, None)
            self.na = NaIA(quaternion, sc_pos, None)
            self.nb = NaIB(quaternion, sc_pos, None)
            self.b0 = BGO0(quaternion, sc_pos, None)
            self.b1 = BGO1(quaternion, sc_pos, None)

        self._detectors = OrderedDict(
            n0=self.n0,
            n1=self.n1,
            n2=self.n2,
            n3=self.n3,
            n4=self.n4,
            n5=self.n5,
            n6=self.n6,
            n7=self.n7,
            n8=self.n8,
            n9=self.n9,
            na=self.na,
            nb=self.nb,
            b0=self.b0,
            b1=self.b1,
        )

    @classmethod
    def from_position_interpolator(
        cls, pos_interp: PositionInterpolator, time: float = 0
    ):
        """
        Create a GBM directly from an interpolator
        """

        met = pos_interp.met(time)

        return cls(pos_interp.quaternion(time), pos_interp.sc_pos(time), gbm_time=met)

    def set_quaternion(self, quaternion):
        """FIXME! briefly describe function

        :param quaternion: 
        :returns: 
        :rtype: 

        """

        for key in self._detectors.keys():
            self._detectors[key].set_quaternion(quaternion)

        self._quaternion = quaternion

    def set_sc_pos(self, sc_pos):
        """FIXME! briefly describe function

        :param sc_pos: 
        :returns: 
        :rtype: 

        """

        for key in self._detectors.keys():
            self._detectors[key].set_sc_pos(sc_pos)

        self._sc_pos = sc_pos

    def get_good_detectors(self, point, fov):
        """
        Returns a list of detectors containing the point in the FOV

        Parameters
        ----------
        point
        fov

        Returns
        -------

        """

        good_detectors = self._contains_point(point, fov)

        return good_detectors

    def get_fov(self, radius, fermi_frame=False):
        """
        Parameters
        ----------
        fermi_frame
        radius

        """

        polys = []

        for key in self._detectors.keys():

            if key[0] == "b":
                this_rad = 90

            else:
                this_rad = radius

            polys.append(self._detectors[key].get_fov(this_rad, fermi_frame))

        polys = np.array(polys)

        return polys

    def get_good_fov(self, point, radius, fermi_frame=False):
        """
        Returns the detectors that contain the given point
        for the given angular radius

        Parameters
        ----------
        point
        radius


        """

        good_detectors = self._contains_point(point, radius)

        polys = []

        for key in good_detectors:
            polys.append(self._detectors[key].get_fov(radius, fermi_frame))

        return [polys, good_detectors]

    def get_sun_angle(self, keys=None):
        """

        Returns
        -------

        """
        angles = []

        if keys is None:

            for key in self._detectors.keys():
                angles.append(self._detectors[key].sun_angle)

        else:

            for key in keys:
                angles.append(self._detectors[key].sun_angle)

        return angles

    def get_centers(self, keys=None):
        """

        Returns
        -------

        """
        centers = []

        if keys is None:

            for key in self._detectors.keys():
                centers.append(self._detectors[key].get_center())

        else:

            for key in keys:
                centers.append(self._detectors[key].get_center())

        return centers

    def get_separation(self, source):
        """
        Get the andular separation of the detectors from a point
        Parameters

        """

        out = {}
        for k, v in self._detectors.items():
            out[k] = v.get_center().separation(source).deg

        return pd.Series(out)

    def get_earth_points(self, fermi_frame=False):
        """

        Returns
        -------

        """

        if self._sc_pos is not None:

            self._calc_earth_points(fermi_frame)

            return self._earth_points

        else:
            print("No spacecraft position set")

    def _calc_earth_points(self, fermi_frame):

        xyz_position = coord.SkyCoord(
            x=self._sc_pos[0],
            y=self._sc_pos[1],
            z=self._sc_pos[2],
            frame="icrs",
            representation="cartesian",
        )

        earth_radius = 6371.0 * u.km
        fermi_radius = np.sqrt((self._sc_pos ** 2).sum())

        horizon_angle = 90 - np.rad2deg(
            np.arccos((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value
        )

        horizon_angle = (180 - horizon_angle) * u.degree

        num_points = 300

        ra_grid_tmp = np.linspace(0, 360, num_points)

        dec_range = [-90, 90]
        cosdec_min = np.cos(np.deg2rad(90.0 + dec_range[0]))
        cosdec_max = np.cos(np.deg2rad(90.0 + dec_range[1]))

        v = np.linspace(cosdec_min, cosdec_max, num_points)
        v = np.arccos(v)
        v = np.rad2deg(v)
        v -= 90.0

        dec_grid_tmp = v

        ra_grid = np.zeros(num_points ** 2)
        dec_grid = np.zeros(num_points ** 2)

        itr = 0
        for ra in ra_grid_tmp:
            for dec in dec_grid_tmp:
                ra_grid[itr] = ra
                dec_grid[itr] = dec
                itr += 1

        if fermi_frame:
            all_sky = coord.SkyCoord(
                Az=ra_grid,
                Zen=dec_grid,
                frame=GBMFrame(quaternion=self._quaternion),
                unit="deg",
            )
        else:
            all_sky = coord.SkyCoord(ra=ra_grid, dec=dec_grid, frame="icrs", unit="deg")

        condition = all_sky.separation(xyz_position) > horizon_angle

        # self.seps = all_sky.separation(xyz_position)

        self._earth_points = all_sky[condition]

    @property
    def detectors(self):

        return self._detectors

    def _contains_point(self, point, radius):
        """
        returns detectors that contain a points
        """

        condition = []

        steps = 500

        for key in self._detectors.keys():

            if key[0] == "b":
                this_rad = 90
            else:
                this_rad = radius

            j2000 = self._detectors[key]._center.icrs

            poly = sp.SphericalPolygon.from_cone(
                j2000.ra.value, j2000.dec.value, this_rad, steps=steps
            )

            if poly.contains_point(point.cartesian.xyz.value):
                condition.append(key)

        return condition

    def plot_detector_pointings(self, ax=None, fov=None, show_earth=True, **kwargs):

        if ax is None:

            skymap_kwargs = {}

            if "projection" in kwargs:
                skymap_kwargs["projection"] = kwargs.pop("projection")

            if "center" in kwargs:
                skymap_kwargs["center"] = kwargs.pop("center")

            if "radius" in kwargs:
                skymap_kwargs["radius"] = kwargs.pop("radius")

            ax = skyplot(**skymap_kwargs)

        _defaults = dict(edgecolor="#13ED9B", lw=1, facecolor="#13ED9B", alpha=0.3)
        for k, v in _defaults.items():
            if k not in kwargs:
                kwargs[k] = v

        for k, v in self._detectors.items():

            v.plot_pointing(ax=ax, fov=fov, **kwargs)

        if show_earth:

            circle = SphericalCircle(
                self._earth_position.icrs.ra,
                self._earth_position.icrs.dec,
                67,
                vertex_unit=u.deg,
                resolution=5000,
                #            edgecolor=color,
                transform=ax.get_transform("icrs"),
                edgecolor="none",
                facecolor="#13ACED",
                alpha=0.1,
                zorder=-3000,
            )

            ax.add_patch(circle)

        return ax.get_figure()

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


# def get_legal_pairs():
#     """
#     Plots the legal pairs of detectors for GBM observations

#     Returns
#     -------

#     """

#     dlp = np.array(
#         [
#             [0, 274, 39, 171, 12, 29, 0, 5, 1, 6, 1, 0],
#             [258, 0, 233, 55, 4, 100, 2, 1, 1, 12, 27, 0],
#             [55, 437, 0, 2, 2, 311, 0, 1, 1, 13, 235, 0],
#             [215, 80, 3, 0, 330, 107, 4, 8, 19, 2, 1, 0],
#             [13, 4, 8, 508, 0, 269, 2, 29, 236, 0, 1, 0],
#             [44, 188, 337, 166, 279, 0, 0, 0, 0, 0, 0, 0],
#             [0, 1, 1, 2, 2, 0, 0, 238, 46, 180, 12, 33],
#             [0, 2, 0, 18, 35, 0, 222, 0, 221, 61, 3, 109],
#             [0, 0, 1, 16, 215, 0, 51, 399, 0, 4, 2, 303],
#             [3, 18, 21, 4, 0, 0, 190, 82, 1, 0, 324, 110],
#             [1, 25, 191, 0, 0, 0, 16, 6, 4, 516, 0, 293],
#             [0, 0, 0, 0, 0, 0, 32, 147, 297, 138, 263, 0],
#         ]
#     )

#     sns.heatmap(dlp, annot=True, fmt="d", cmap="YlGnBu")
#     plt.ylabel("NaI")
#     plt.xlabel("NaI")
