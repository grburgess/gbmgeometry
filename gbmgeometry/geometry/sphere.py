import ipyvolume as ipv
import numpy as np
import h5py

from gbmgeometry.utils.package_utils import get_path_of_data_file


# -*- coding: utf-8 -*-
class Sphere(object):
    def __init__(
        self, ax, x=0, y=0, z=0, radius=1.0, detail_level=100, color="#FFFFFF", **kwargs
    ):

        self._x = x
        self._y = y
        self._z = z

        self._radius = radius

        self._detail_level = detail_level
        self._color = color

    @property
    def radius(self):
        return self._radius

    def plot(self, **kwargs):
        """
        
        plot the sphere
        

        :returns: 
        :rtype: 

        """

        u = np.linspace(0, 2 * np.pi, self._detail_level)
        v = np.linspace(0, np.pi, self._detail_level)

        if np.atleast_1d(self._x).shape[0] == 1:

            X = self._x + self._radius * np.outer(np.cos(u), np.sin(v))

            Y = self._y + self._radius * np.outer(np.sin(u), np.sin(v))

            Z = self._z + self._radius * np.outer(np.ones(np.size(u)), np.cos(v))

        else:

            # for animations

            X = np.array(
                [x + self._radius * np.outer(np.cos(u), np.sin(v)) for x in self._x]
            )

            Y = np.array(
                [y + self._radius * np.outer(np.sin(u), np.sin(v)) for y in self._y]
            )

            Z = np.array(
                [
                    z + self._radius * np.outer(np.ones(np.size(u)), np.cos(v))
                    for z in self._z
                ]
            )

        return ipv.plot_surface(X, Y, Z, color=self._color, **kwargs)
