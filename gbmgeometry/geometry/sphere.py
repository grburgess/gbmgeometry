import ipyvolume as ipv
import numpy as np
import matplotlib.pyplot as plt

import h5py

from gbmgeometry.utils.package_utils import get_path_of_data_file


# -*- coding: utf-8 -*-
class Sphere(object):
    def __init__(
        self,
        ax,
        x=0,
        y=0,
        z=0,
        radius=1.0,
        detail_level=100,
        color="#FFFFFF",
        image=None,
        **kwargs
    ):

        self._x = x
        self._y = y
        self._z = z

        self._radius = radius

        self._detail_level = detail_level
        self._color = color

        self._image = image

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

        x_unit = np.outer(np.cos(u), np.sin(v))
        y_unit = np.outer(np.sin(u), np.sin(v))
        z_unit = np.outer(np.ones(np.size(u)), np.cos(v))

        if np.atleast_1d(self._x).shape[0] == 1:

            X = self._x + self._radius * x_unit

            Y = self._y + self._radius * y_unit

            Z = self._z + self._radius * z_unit

        else:

            X = np.array([x + self._radius * x_unit for x in self._x])

            Y = np.array([y + self._radius * y_unit for y in self._y])

            Z = np.array([z + self._radius * z_unit for z in self._z])

        if self._image is None:

            return ipv.plot_surface(X, Y, Z, color=self._color, **kwargs)

        else:

            lon = np.arctan2(y_unit, x_unit)
            lat = np.arcsin(z_unit)

            u = 0.5 + lon / (2 * np.pi)
            v = 0.5 + lat / (np.pi)

            return ipv.plot_mesh(
                X, Y, Z, u=u, v=v, texture=self._image, wireframe=False
            )
