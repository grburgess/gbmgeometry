from .gbm_detector import NaI0, NaI1, NaI2, NaI3, NaI4, NaI5
from .gbm_detector import NaI6, NaI7, NaI8, NaI9, NaIA, NaIB
from .gbm_detector import BGO0, BGO1
import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt

import numpy as np
from collections import OrderedDict
import spherical_geometry.polygon as sp

from astropy.table import Table
import astropy.units as u

import seaborn as sns

_det_color_cycle = np.linspace(0, 1, 12)


class GBM(object):
    def __init__(self, quaternion, sc_pos=None):

        """

        Parameters
        ----------
        quaternion : Fermi GBM quarternion array
        """
        self.n0 = NaI0(quaternion, sc_pos)
        self.n1 = NaI1(quaternion, sc_pos)
        self.n2 = NaI2(quaternion, sc_pos)
        self.n3 = NaI3(quaternion, sc_pos)
        self.n4 = NaI4(quaternion, sc_pos)
        self.n5 = NaI5(quaternion, sc_pos)
        self.n6 = NaI6(quaternion, sc_pos)
        self.n7 = NaI7(quaternion, sc_pos)
        self.n8 = NaI8(quaternion, sc_pos)
        self.n9 = NaI9(quaternion, sc_pos)
        self.na = NaIA(quaternion, sc_pos)
        self.nb = NaIB(quaternion, sc_pos)
        self.b0 = BGO0(quaternion, sc_pos)
        self.b1 = BGO1(quaternion, sc_pos)

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

    def set_quaternion(self, quaternion):
        """
        Parameters
        ----------
        quaternion
        """
        for key in self._detectors.keys():
            self._detectors[key].set_quaternion(quaternion)

    def set_sc_pos(self, sc_pos):
        """
        Parameters
        ----------
        sc_pos
        """

        for key in self._detectors.keys():
            self._detectors[key].set_sc_pos(sc_pos)

    def get_fov(self, radius):
        """
        Parameters
        ----------
        radius

        """

        polys = []

        for key in self._detectors.keys():

            if key[0] == 'b':
                this_rad = 120

            else:
                this_rad = radius

            polys.append(self._detectors[key].get_fov(this_rad))

        polys = np.array(polys)

        return polys

    def get_good_fov(self, point, radius):
        """
        Returns the detectors that contain the given point
        for the given angular radius

        Parameters
        ----------
        point
        radius


        """

        good_detectors = self._contains_point(point, radius)

        polys = self.get_fov(radius)

        return [polys[good_detectors], np.where(good_detectors)[0]]

    def get_centers(self):

        """

        Returns
        -------

        """
        centers = []
        for key in self._detectors.keys():
            centers.append(self._detectors[key].get_center())

        return centers

    def detector_plot(self, radius=60., point=None, good=False, projection='moll', lat_0=0, lon_0=0, fignum=1,
                      map=None):

        """

        Parameters
        ----------
        radius
        point
        good
        projection
        lat_0
        lon_0
        """

        map_flag = False
        if map is None:

            fig = plt.figure(fignum)
            ax = fig.add_subplot(111)

            map = bm.Basemap(projection=projection, lat_0=lat_0, lon_0=lon_0,
                             resolution='l', area_thresh=1000.0, celestial=True, ax=ax)


        else:

            map_flag = True

        good_detectors = range(14)
        centers = self.get_centers()

        if good and point:

            fovs, good_detectors = self.get_good_fov(point, radius)
            map.plot(point.ra.value, point.dec.value, '*', color='#ffffbf', latlon=True)




        else:

            fovs = self.get_fov(radius)

        if point:
            map.plot(point.ra.value, point.dec.value, '*', color='#ffffbf', latlon=True)

        color_itr = np.linspace(0, 1, len(fovs))

        for i, fov in enumerate(fovs):
            ra, dec = fov

            idx = np.argsort(ra)

            map.plot(ra[idx], dec[idx], '.', color=plt.cm.Set1(color_itr[i]), latlon=True, markersize=2.)

            x, y = map(centers[good_detectors[i]].icrs.ra.value, centers[good_detectors[i]].icrs.dec.value)

            plt.text(x, y, self._detectors.keys()[good_detectors[i]], color=plt.cm.Set1(color_itr[i]))

        if not map_flag:
            _ = map.drawmeridians(np.arange(0, 360, 30), color='#3A3A3A')
            _ = map.drawparallels(np.arange(-90, 90, 15), labels=[True] * len(np.arange(-90, 90, 15)), color='#3A3A3A')
            map.drawmapboundary(fill_color='#151719')

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

    def _contains_point(self, point, radius):
        """
        returns detectors that contain a points
        """

        condition = []

        steps = 500

        for key in self._detectors.keys():

            if key[0] =='b':
                this_rad = 180
            else:
                this_rad = radius

            j2000 = self._detectors[key]._center.icrs

            poly = sp.SphericalPolygon.from_cone(j2000.ra.value,
                                                 j2000.dec.value,
                                                 this_rad,
                                                 steps=steps)

            condition.append(poly.contains_point(point.cartesian.xyz.value))

        return np.array(condition)


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
