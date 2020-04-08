import ipyvolume as ipv
import numpy as np
import matplotlib.pyplot as plt
import numba as nb
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
        transform_matrix=None,
        **kwargs,
    ):

        self._x = x
        self._y = y
        self._z = z

        self._radius = radius

        self._detail_level = int(detail_level)
        self._color = color

        self._image = image
        self._transform_matrix = transform_matrix

        self._time_dep_transform = False

        if (transform_matrix is not None) and (len(transform_matrix.shape) == 3):
            self._time_dep_transform = True
            self._n_steps = transform_matrix.shape[0]

    @property
    def radius(self):
        return self._radius

    def plot(self, **kwargs):
        """
        
        plot the sphere
        

        :returns: 
        :rtype: 

        """

        # u = np.linspace(0, 2 * np.pi, self._detail_level)
        # v = np.linspace(0, np.pi, self._detail_level)

        # x_unit = np.outer(np.cos(u), np.sin(v))
        # y_unit = np.outer(np.sin(u), np.sin(v))
        # z_unit = np.outer(np.ones(np.size(u)), np.cos(v))

        u = np.linspace(0, 1, self._detail_level)
        v = np.linspace(0, 1, self._detail_level)
        u, v = np.meshgrid(u, v)
        phi = u * 2 * np.pi
        theta = v * np.pi


        x_unit = np.cos(phi) * np.sin(theta)
        y_unit = np.sin(theta) * np.sin(phi)
        z_unit = np.cos(theta)

        
        
        if self._transform_matrix is not None:

            xyz = np.array([x_unit, y_unit, z_unit]).T

            if self._time_dep_transform:

                # new_xyz = compute_multiple_rotation(xyz, self._transform_matrix, self._detail_level, self._n_steps)

                x_unit_rot = np.zeros(
                    (self._n_steps, self._detail_level, self._detail_level)
                )
                y_unit_rot = np.zeros(
                    (self._n_steps, self._detail_level, self._detail_level)
                )
                z_unit_rot = np.zeros(
                    (self._n_steps, self._detail_level, self._detail_level)
                )

                for i in range(self._n_steps):

                    this_xyz = compute_single_rotation(
                        xyz, self._transform_matrix[i], self._detail_level
                    )

                    x_unit_rot[i] = this_xyz[0]
                    y_unit_rot[i] = this_xyz[1]
                    z_unit_rot[i] = this_xyz[2]

            else:

                xyz = compute_single_rotation(
                    xyz, self._transform_matrix, self._detail_level
                )

                x_unit_rot = xyz[0]
                y_unit_rot = xyz[1]
                z_unit_rot = xyz[2]

        if np.atleast_1d(self._x).shape[0] == 1:

            if self._time_dep_transform:
                # if False:
                X = np.array(
                    [self._x + self._radius * x_unit for _ in range(self._n_steps)]
                )

                Y = np.array(
                    [self._y + self._radius * y_unit for _ in range(self._n_steps)]
                )

                Z = np.array(
                    [self._z + self._radius * z_unit for _ in range(self._n_steps)]
                )

            else:

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

            if self._transform_matrix is None:

                lon = np.arctan2(y_unit, x_unit)
                lat = np.arcsin(z_unit)

            else:

                lon = np.arctan2(y_unit_rot, x_unit_rot)
                lat = np.arcsin(z_unit_rot)

            u = 0.5 + lon / (2 * np.pi)
            v = 0.5 + lat / (np.pi)

            return ipv.plot_mesh(
                X, Y, Z, u=u, v=v, texture=self._image, wireframe=False
            )


@nb.njit(fastmath=True)
def compute_single_rotation(xyz, transform_matrix, detail_level):

    new_xyz = np.zeros((detail_level, detail_level, 3))

    for i in range(detail_level):
        for j in range(detail_level):
            new_xyz[i, j] = np.dot(transform_matrix, xyz[i, j, :])

    return new_xyz.T


@nb.njit(fastmath=True)
def compute_multiple_rotation(xyz, transform_matrix, detail_level, time_steps):

    new_xyz = np.zeros((time_steps, detail_level, detail_level, 3))
    #    out_xyz = np.zeros((time_steps, 3, detail_level, detail_level))

    for i in range(time_steps):
        for j in range(detail_level):
            for k in range(detail_level):
                new_xyz[i, j, k] = np.dot(transform_matrix[i], xyz[j, k, :])

    #        out_xyz[i] = new_xyz[i].T

    return new_xyz
