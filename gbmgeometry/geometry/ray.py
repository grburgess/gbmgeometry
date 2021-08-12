import collections
import numpy as np
from sympy import Point3D, Line3D


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
        theta = np.deg2rad(90.0 - self._point_source.spherical.lat.value)
        phi = np.deg2rad(self._point_source.spherical.lon.value)

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
