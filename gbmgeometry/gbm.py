from .gbm_detector import NaI0, NaI1, NaI2, NaI3, NaI4, NaI5
from .gbm_detector import NaI6, NaI7, NaI8, NaI9, NaIA, NaIB

import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt

import numpy as np
from collections import OrderedDict
from spherical_geometry.polygon import SphericalPolygon

_det_color_cycle = np.linspace(0, 1, 12)


class GBM(object):
    def __init__(self, quaternion):

        """

        Parameters
        ----------
        quaternion : Fermi GBM quarternion array
        """
        self.n0 = NaI0(quaternion)
        self.n1 = NaI1(quaternion)
        self.n2 = NaI2(quaternion)
        self.n3 = NaI3(quaternion)
        self.n4 = NaI4(quaternion)
        self.n5 = NaI5(quaternion)
        self.n6 = NaI6(quaternion)
        self.n7 = NaI7(quaternion)
        self.n8 = NaI8(quaternion)
        self.n9 = NaI9(quaternion)
        self.na = NaIA(quaternion)
        self.nb = NaIB(quaternion)

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
                                      nb=self.nb)

    def set_quarternion(self, quaternion):
        """
        """
        for key in self._detectors.keys():
            self._detectors[key].set_quarternion(quaternion)

    def get_fov(self, radius):
        """

        """

        polys = []

        for key in self._detectors.keys():
            polys.append(self._detectors[key].get_fov(radius))

        polys = np.array(polys)

        return polys

    def get_good_fov(self, point, radius):
        """
        Returns the detectors that contain the given point
        for the given angular radius


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

    def detector_plot(self, radius=60., point=None, good=False, projection='moll', lat_0=0, lon_0=0):

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
        map = bm.Basemap(projection=projection, lat_0=lat_0, lon_0=lon_0,
                         resolution='l', area_thresh=1000.0, celestial=True)

        map.drawmapboundary(fill_color='#5B5655')

        good_detectors = range(12)
        centers = self.get_centers()

        if good and point:

            fovs, good_detectors = self.get_good_fov(point, radius)
            map.plot(point.ra.value, point.dec.value, '*', color='yellow', latlon=True)




        else:

            fovs = self.get_fov(radius)

        if point:
            map.plot(point.ra.value, point.dec.value, '*', color='yellow', latlon=True)

        color_itr = np.linspace(0, 1, len(fovs))

        for i, fov in enumerate(fovs):
            ra, dec = fov

            map.plot(ra, dec, '.', color=plt.cm.Set1(color_itr[i]), latlon=True, markersize=2.)

            x, y = map(centers[good_detectors[i]].icrs.ra.value, centers[good_detectors[i]].icrs.dec.value)

            plt.text(x, y, self._detectors.keys()[good_detectors[i]], color='w')

        _ = map.drawmeridians(np.arange(0, 360, 30), color='#3A3A3A')
        _ = map.drawparallels(np.arange(-90, 90, 15), labels=[True] * len(np.arange(-90, 90, 15)), color='#3A3A3A')

    def _contains_point(self, point, radius):
        """
        returns detectors that contain a points
        """

        condition = []

        steps = 500

        for key in self._detectors.keys():
            j2000 = self._detectors[key]._center.icrs

            poly = SphericalPolygon.from_cone(j2000.ra.value,
                                              j2000.dec.value,
                                              radius,
                                              steps=steps)

            condition.append(poly.contains_point(point.cartesian.xyz.value))

        return np.array(condition)
