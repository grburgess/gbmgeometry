from collections import OrderedDict

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
#import mpl_toolkits.basemap as bm
import numpy as np
import spherical_geometry.polygon as sp
from astropy.table import Table
import astropy.time as time

from .gbm_detector import BGO0, BGO1
from .gbm_detector import NaI0, NaI1, NaI2, NaI3, NaI4, NaI5
from .gbm_detector import NaI6, NaI7, NaI8, NaI9, NaIA, NaIB
from .gbm_frame import GBMFrame


from gbmgeometry.utils.gbm_time import GBMTime

# import seaborn as sns

_det_color_cycle = np.linspace(0, 1, 12)


class GBM(object):
    def __init__(self, quaternion, sc_pos=None, gbm_time=None):

        """

        Parameters
        ----------
        quaternion : Fermi GBM quarternion array
        """

        if gbm_time is not None:

            if isinstance(gbm_time, str):

               self._gbm_time = GBMTime.from_UTC_fits(gbm_time)

            else:

                # assuming MET

                self._gbm_time = GBMTime.from_MET(gbm_time)

        else:

            self._gbm_time = None







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

        self._detectors = OrderedDict(n0=self.n0,
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
                                      b1=self.b1)
        self._quaternion = quaternion
        self._sc_pos = sc_pos

    def set_quaternion(self, quaternion):
        """
        Parameters
        ----------
        quaternion
        """
        for key in self._detectors.keys():
            self._detectors[key].set_quaternion(quaternion)

        self._quaternion = quaternion

    def set_sc_pos(self, sc_pos):
        """
        Parameters
        ----------
        sc_pos
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

            if key[0] == 'b':
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

    # def plot_pointing(self, ra_0=0, dec_0=0, projection='moll', fignum=1, point=None):
    #     """
    #
    #
    #     Returns
    #     -------
    #
    #     """
    #
    #     fig = plt.figure(fignum)
    #     ax = fig.add_subplot(111)
    #
    #     map = bm.Basemap(projection=projection, lat_0=dec_0, lon_0=ra_0, celestial=False, ax=ax)
    #
    #     color_itr = np.linspace(0, .8, 14)
    #
    #     for i, center in enumerate(self.get_centers()):
    #         ra, dec = center.icrs.ra.value, center.icrs.dec.value
    #
    #         idx = np.argsort(ra)
    #
    #         map.plot(ra, dec, '.', color=plt.cm.Set2(color_itr[i]), latlon=True, markersize=3.)
    #
    #     if point is not None:
    #         ra, dec = point.icrs.ra.value, point.icrs.dec.value
    #
    #         map.plot(ra, dec, '*', color='yellow', latlon=True, markersize=3.)
    #
    #     _ = map.drawmeridians(np.arange(0, 360, 30), color='#3A3A3A')
    #     _ = map.drawparallels(np.arange(-90, 90, 15), labels=[True] * len(np.arange(-90, 90, 15)), color='#3A3A3A')
    #     map.drawmapboundary(fill_color='#151719')

    # def detector_plot(self, radius=60., point=None, good=False, projection='moll', lat_0=0, lon_0=0, fignum=1,
    #                   map=None, show_earth=False, fermi_frame=False):
    #
    #     """
    #
    #     Parameters
    #     ----------
    #     radius
    #     point
    #     good
    #     projection
    #     lat_0
    #     lon_0
    #     """
    #
    #     map_flag = False
    #     if map is None:
    #
    #         fig = plt.figure(fignum)
    #         ax = fig.add_subplot(111)
    #
    #         map = bm.Basemap(projection=projection, lat_0=lat_0, lon_0=lon_0,
    #                          resolution='l', area_thresh=1000.0, celestial=True, ax=ax)
    #
    #
    #     else:
    #
    #         map_flag = True
    #
    #     good_detectors = self._detectors.keys()
    #
    #     if good and point:
    #
    #         fovs, good_detectors = self.get_good_fov(point, radius, fermi_frame)
    #         # map.plot(point.ra.value, point.dec.value, '*', color='#ffffbf', latlon=True)
    #
    #
    #
    #
    #     else:
    #
    #         fovs = self.get_fov(radius, fermi_frame)
    #
    #     if point:
    #         pass
    #         # map.plot(point.ra.value, point.dec.value, '*', color='#ffffbf', latlon=True)
    #
    #     color_itr = np.linspace(0, .8, len(fovs))
    #
    #     for i, fov in enumerate(fovs):
    #         ra, dec = fov
    #
    #         idx = np.argsort(ra)
    #
    #         map.plot(ra[idx], dec[idx], '.', color=plt.cm.Set2(color_itr[i]), latlon=True, markersize=4.)
    #
    #         if fermi_frame:
    #             x, y = map(self._detectors[good_detectors[i]].get_center().Az.value,
    #                        self._detectors[good_detectors[i]].get_center().Zen.value)
    #         else:
    #             x, y = map(self._detectors[good_detectors[i]].get_center().icrs.ra.value,
    #                        self._detectors[good_detectors[i]].get_center().icrs.dec.value)
    #
    #         plt.text(x, y, good_detectors[i], color=plt.cm.Set2(color_itr[i]), size=9)
    #
    #     if show_earth and self._sc_pos is not None:
    #
    #         earth_points = self.get_earth_points(fermi_frame)
    #
    #         if fermi_frame:
    #
    #             lon, lat = earth_points.Az.value, earth_points.Zen.value
    #
    #         else:
    #
    #             lon, lat = earth_points.ra.value, earth_points.dec.value
    #
    #         idx = np.argsort(lon)
    #         lon = lon[idx]
    #         lat = lat[idx]
    #
    #         map.plot(lon, lat, '.', color="#0C81F9", latlon=True, alpha=0.35, markersize=4.5)
    #
    #     if not map_flag:
    #         _ = map.drawmeridians(np.arange(0, 360, 30), color='#3A3A3A')
    #         _ = map.drawparallels(np.arange(-90, 90, 15), labels=[True] * len(np.arange(-90, 90, 15)), color='#3A3A3A')
    #         map.drawmapboundary(fill_color='#151719')

    def get_separation(self, source):
        """
        Get the andular separation of the detectors from a point
        Parameters
        ----------
        source

        Returns
        -------

        """

        tab = Table(names=["Detector", "Separation"], dtype=["|S2", np.float64])

        for key in self._detectors.keys():
            sep = self._detectors[key].get_center().separation(source)
            tab.add_row([key, sep])

        tab['Separation'].unit = u.degree

        tab.sort("Separation")

        return tab

    def get_earth_points(self, fermi_frame=False):
        """

        Returns
        -------

        """

        if self._sc_pos is not None:

            self._calc_earth_points(fermi_frame)

            return self._earth_points

        else:
            print "No spacecraft position set"

    def _calc_earth_points(self, fermi_frame):

        xyz_position = coord.SkyCoord(x=self._sc_pos[0],
                                      y=self._sc_pos[1],
                                      z=self._sc_pos[2],
                                      frame='icrs',
                                      representation='cartesian')

        earth_radius = 6371. * u.km
        fermi_radius = np.sqrt((self._sc_pos ** 2).sum())

        horizon_angle = 90 - np.rad2deg(np.arccos((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value)

        horizon_angle = (180 - horizon_angle) * u.degree

        num_points = 300

        ra_grid_tmp = np.linspace(0, 360, num_points)

        dec_range = [-90, 90]
        cosdec_min = np.cos(np.deg2rad(90.0 + dec_range[0]))
        cosdec_max = np.cos(np.deg2rad(90.0 + dec_range[1]))

        v = np.linspace(cosdec_min, cosdec_max, num_points)
        v = np.arccos(v)
        v = np.rad2deg(v)
        v -= 90.

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
            all_sky = coord.SkyCoord(Az=ra_grid, Zen=dec_grid, frame=GBMFrame(quaternion=self._quaternion), unit='deg')
        else:
            all_sky = coord.SkyCoord(ra=ra_grid, dec=dec_grid, frame='icrs', unit='deg')

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

            if key[0] == 'b':
                this_rad = 90
            else:
                this_rad = radius

            j2000 = self._detectors[key]._center.icrs

            poly = sp.SphericalPolygon.from_cone(j2000.ra.value,
                                                 j2000.dec.value,
                                                 this_rad,
                                                 steps=steps)

            if poly.contains_point(point.cartesian.xyz.value):
                condition.append(key)

        return condition


def get_legal_pairs():
    """
    Plots the legal pairs of detectors for GBM observations

    Returns
    -------

    """

    dlp = np.array([[0, 274, 39, 171, 12, 29, 0, 5, 1, 6, 1, 0],
                    [258, 0, 233, 55, 4, 100, 2, 1, 1, 12, 27, 0],
                    [55, 437, 0, 2, 2, 311, 0, 1, 1, 13, 235, 0],
                    [215, 80, 3, 0, 330, 107, 4, 8, 19, 2, 1, 0],
                    [13, 4, 8, 508, 0, 269, 2, 29, 236, 0, 1, 0],
                    [44, 188, 337, 166, 279, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 2, 2, 0, 0, 238, 46, 180, 12, 33],
                    [0, 2, 0, 18, 35, 0, 222, 0, 221, 61, 3, 109],
                    [0, 0, 1, 16, 215, 0, 51, 399, 0, 4, 2, 303],
                    [3, 18, 21, 4, 0, 0, 190, 82, 1, 0, 324, 110],
                    [1, 25, 191, 0, 0, 0, 16, 6, 4, 516, 0, 293],
                    [0, 0, 0, 0, 0, 0, 32, 147, 297, 138, 263, 0]])

    sns.heatmap(dlp, annot=True, fmt='d', cmap="YlGnBu")
    plt.ylabel("NaI")
    plt.xlabel("NaI")
