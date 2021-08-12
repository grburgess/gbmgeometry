import collections

import ipyvolume as ipv
import numpy as np
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sympy import Line3D, Plane, Point3D

from gbmgeometry.geometry.cube import Cube


def get_sc_matrix(q1, q2, q3, q4):
    """
    build a sc matrix from quats

    :param q1: 
    :param q2: 
    :param q3: 
    :param q4: 
    :returns: 
    :rtype: 

    """

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

    return sc_matrix


class Ray(object):
    # scaling for distance of lines
    _R = 1e10  # cm
    _scale = 3e-8

    def __init__(self, detector, point_source, color="#29FC5C", probability=None):
        self._detector = detector

        self._probability = probability

        self._color = color

        self._point_source = point_source

        self._calculate_ray_origin()

    def _calculate_ray_origin(self):
        theta = np.deg2rad(90.0 - self._point_source.lat.value)
        phi = np.deg2rad(self._point_source.lon.value)

        x = Ray._R * np.cos(phi) * np.sin(theta)
        y = Ray._R * np.sin(phi) * np.sin(theta)
        z = Ray._R * np.cos(theta)

        # this is the "distant orgin of the ray"
        self._origin = np.array([x, y, z])

        self._sympy_line = Line3D(
            Point3D(self._detector.mount_point), Point3D(self._origin)
        )

        self._plot_origin = self.point_on_ray(Ray._scale)

    def plot(self, ax):
        ax.plot(
            [self._plot_origin[0], self.detector_origin[0]],
            [self._plot_origin[1], self.detector_origin[1]],
            [self._plot_origin[2], self.detector_origin[2]],
            color=self._color,
            alpha=0.8,
        )

    @property
    def detector_name(self):
        return self._detector.name

    @property
    def probability(self):
        return self._probability

    @property
    def detector_origin(self):
        return self._detector.mount_point

    @property
    def ray_origin(self):
        return self._origin

    @property
    def sympy_line(self):
        return self._sympy_line

    def point_on_ray(self, t=0.5):
        """
        get a parametrized point on the ray
        
        :param t: point between 0 and 1
        :return: 
        """

        assert 0.0 <= t <= 1.0, "t must be between 0 and 1"

        return self.detector_origin + (self._origin - self.detector_origin) * t


class Surface(object):
    def __init__(self, name, vertices):

        self._vertices = vertices

        self._name = name

        if "+" in name:

            self._sign = 1.0

        elif "-" in name:

            self._sign = -1.0

        else:

            raise RuntimeError("the plane name is wrong")

        self._calculate_origin()

    def _calculate_origin(self):
        """
        
        compute the origin of the plane
        
        
        :return: 
        """

        if "x" in self._name:

            assert len(np.unique(self._vertices[:, 0])) == 1, "vertices are wrong!"

            self._normal = np.array([1, 0, 0]) * self._sign

            x_origin = self._vertices[0, 0]

            y_origin = (min(self._vertices[:, 1]) + max(self._vertices[:, 1])) / 2.0

            z_origin = (min(self._vertices[:, 2]) + max(self._vertices[:, 2])) / 2.0

            # self._edges = [self.]

        elif "y" in self._name:

            assert len(np.unique(self._vertices[:, 1])) == 1, "vertices are wrong!"

            self._normal = np.array([0, 1.0, 0]) * self._sign

            x_origin = (min(self._vertices[:, 0]) + max(self._vertices[:, 0])) / 2.0

            y_origin = self._vertices[0, 1]

            z_origin = (min(self._vertices[:, 2]) + max(self._vertices[:, 2])) / 2.0

        elif "z" in self._name:

            assert len(np.unique(self._vertices[:, 2])) == 1, "vertices are wrong!"

            self._normal = np.array([0, 0, 1]) * self._sign

            x_origin = (min(self._vertices[:, 0]) + max(self._vertices[:, 0])) / 2.0

            y_origin = (min(self._vertices[:, 1]) + max(self._vertices[:, 1])) / 2.0

            z_origin = self._vertices[0, 2]

        self._xmax = self._vertices[:, 0].max()
        self._ymax = self._vertices[:, 1].max()
        self._zmax = self._vertices[:, 2].max()

        self._xmin = self._vertices[:, 0].min()
        self._ymin = self._vertices[:, 1].min()
        self._zmin = self._vertices[:, 2].min()

        self._origin = np.array([x_origin, y_origin, z_origin])

        # TODO: perhaps rewrite this to find the plane

        self._sympy_plane = Plane(Point3D(self._origin), normal_vector=self._normal)

    def is_intersecting(self, ray):
        """
        checks if ray intersects plane
        
        :param ray: 
        :return: bool, array
        """

        intersecting_point = self._sympy_plane.intersection(ray.sympy_line)[0]

        if "x" in self._name:

            if self._within_y_bounds(intersecting_point.y) and self._within_z_bounds(
                intersecting_point.z
            ):
                return (
                    True,
                    np.array(
                        [float(x) for x in
                            [
                                intersecting_point.x,
                                intersecting_point.y,
                                intersecting_point.z,
                            ]
                        ]
                    ),
                )

        elif "y" in self._name:

            if self._within_x_bounds(intersecting_point.x) and self._within_z_bounds(
                intersecting_point.z
            ):
                return (
                    True,
                    np.array(
                        [float(x) for x in
                            [
                                intersecting_point.x,
                                intersecting_point.y,
                                intersecting_point.z,
                            ]
                        ]
                    ),
                )

        elif "z" in self._name:

            if self._within_y_bounds(intersecting_point.y) and self._within_x_bounds(
                intersecting_point.x
            ):
                return (
                    True,
                    np.array(
                          [float(x) for x in
                            [
                                intersecting_point.x,
                                intersecting_point.y,
                                intersecting_point.z,
                            ]
                           ]
                    ),
                )

        return False, None

    def _within_x_bounds(self, x):

        if x <= self._xmax and self._xmin <= x:

            return True

        else:

            return False

    def _within_y_bounds(self, y):

        if y <= self._ymax and self._ymin <= y:

            return True

        else:

            return False

    def _within_z_bounds(self, z):

        if z <= self._zmax and self._zmin <= z:

            return True

        else:

            return False

    @property
    def origin(self):

        return self._origin


class Volume(object):
    def __init__(
        self,
        name,
        x_origin,
        y_origin,
        z_origin,
        height,
        x_width,
        y_width,
        transform_matrix=None,
        sc_pos=None,
        color="grey",
        active_surfaces=None,
    ):

        if active_surfaces is None:

            self._active_surfaces = ("+x", "-x", "+y", "-y", "+z", "-z")

        else:

            for surface in active_surfaces:
                assert surface in [
                    "+x",
                    "-x",
                    "+y",
                    "-y",
                    "+z",
                    "-z",
                ], "not a valid surface!"

            self._active_surfaces = active_surfaces

        self._center = (x_origin, y_origin, z_origin)

        self._color = colors.to_rgb(color)

        self._name = name

        # this is to transform into space coords

        # if transform_matrix is not None:

        #     x_origin = x_origin + sc_pos[0]
        #     y_origin = y_origin + sc_pos[1]
        #     z_origin = z_origin + sc_pos[2]

        self._sc_pos = sc_pos
        self._transform_matrix = transform_matrix

        self._x_origin = x_origin - x_width / 2.0
        self._y_origin = z_origin - height / 2.0
        self._z_origin = y_origin - y_width / 2.0

        self._build_cube(
            origin=(
                x_origin - x_width / 2.0,
                z_origin - height / 2.0,
                y_origin - y_width / 2.0,
            ),
            width=x_width,
            depth=y_width,
            height=height,
        )

        self._x_width = x_width
        self._y_width = y_width
        self._z_width = height
        self._intersections = None

    def _build_cube(self, origin=None, width=1, height=1, depth=1):

        self._planes = collections.OrderedDict()
        for plane in ["+x", "-x", "+y", "-y", "+z", "-z"]:
            self._planes[plane] = None

        u, v, w = (0, 0, 0) if origin is None else origin

        grids = []

        for plane in self._planes.keys():

            if "-z" in plane:
                this_grid = self._grid("xy", (u, w), width, depth, v)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

            if "+z" in plane:
                this_grid = self._grid("xy", (u, w), width, depth, v + height)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

            if "-y" in plane:
                this_grid = self._grid("xz", (u, v), width, height, w)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

            if "+y" in plane:
                this_grid = self._grid("xz", (u, v), width, height, w + depth)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

            if "-x" in plane:
                this_grid = self._grid("yz", (w, v), depth, height, u)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

            if "+x" in plane:
                this_grid = self._grid("yz", (w, v), depth, height, u + width)

                self._planes[plane] = Surface(plane, this_grid[0])

                grids.extend(this_grid)

        self._quads = np.array(grids)

    @staticmethod
    def _grid(plane="xy", origin=None, width=1, height=1, depth=0):

        u, v = (0, 0) if origin is None else origin

        # hard code this
        width_segments = 1
        height_segments = 1

        w_x, h_y = width, height

        quads = []

        for i in range(width_segments):
            for j in range(height_segments):
                quads.append(
                    Volume._quad(plane, (i * w_x + u, j * h_y + v), w_x, h_y, depth)
                )

        return np.array(quads)

    @staticmethod
    def _quad(plane="xy", origin=None, width=1, height=1, depth=0):
        u, v = (0, 0) if origin is None else origin

        plane = plane.lower()
        if plane == "xy":
            vertices = (
                (u, v, depth),
                (u + width, v, depth),
                (u + width, v + height, depth),
                (u, v + height, depth),
            )
        elif plane == "xz":
            vertices = (
                (u, depth, v),
                (u + width, depth, v),
                (u + width, depth, v + height),
                (u, depth, v + height),
            )
        elif plane == "yz":
            vertices = (
                (depth, u, v),
                (depth, u + width, v),
                (depth, u + width, v + height),
                (depth, u, v + height),
            )
        else:
            raise ValueError('"{0}" is not a supported plane!'.format(plane))

        return np.array(vertices)

    @property
    def planes(self):
        return self._planes

    @property
    def center(self):
        return self._center

    @property
    def name(self):

        return self._name

    def plot_ipv(self):

        x, y, z = self._center

        cube = Cube(
            color=self._color,
            # x=self._x_origin,
            # y=self._z_origin,
            # z=self._y_origin,
            x=x,
            y=y,
            z=z,
            x_width=self._x_width,
            z_width=self._z_width,
            y_width=self._y_width,
            transform_matrix=self._transform_matrix,
            sc_pos=self._sc_pos,
        )
        return cube.plot()

    def plot(self, ax, alpha=0.1):

        collection = Poly3DCollection(self._quads, facecolors=self._color, alpha=0.25)

        c = []

        for i in self._color:
            c.append(i)

        c.append(alpha)

        collection.set_facecolor(c)

        ax.add_collection3d(collection)

    def intersect_ray(self, ray):

        intersections = collections.OrderedDict()

        for k, v in self._planes.items():

            if k in self._active_surfaces:

                intersection_info = collections.OrderedDict()

                is_intersecting, point = v.is_intersecting(ray)

                if is_intersecting:
                    intersection_info["intersection point"] = point

                    # now get the distance between the points

                    d2 = (np.power(point - ray.detector_origin, 2)).sum()

                    intersection_info["distance"] = np.sqrt(d2)

                    intersections[k] = intersection_info

        self._intersections = intersections

    @property
    def all_intersections(self):

        return self._intersections

    @property
    def intersection(self):

        # return the first intersection

        max_distance = 0.0

        intersection = None

        for k, v in self._intersections.items():

            if v["distance"] > max_distance:
                intersection = k
                max_distance = v["distance"]

        if intersection is None:
            return None, None, None

        return (
            intersection,
            self._intersections[intersection]["intersection point"],
            self._intersections[intersection]["distance"],
        )


# -*- coding: utf-8 -*-
