import numpy as np
from sympy import Plane, Point3D, Line3D


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
