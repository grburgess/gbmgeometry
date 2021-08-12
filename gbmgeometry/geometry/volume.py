import numpy as np
import collections
from gbmgeometry.geometry import Surface
from gbmgeometry.geometry import Cube
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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
        quaternion=None,
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
        self._quaternion = quaternion

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
            quarternion=self._quaternion,
        )

        cube.plot()

        return cube.artists

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
