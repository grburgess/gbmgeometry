import numpy as np
import ipyvolume as ipv


class Cone(object):
    def __init__(self, xo, yo, zo, xp, yp, zp, angle, detail_level=20):

        self._xo = xo
        self._yo = yo
        self._zo = zo

        self._xp = xp
        self._yp = yp
        self._zp = zp

        self._angle = np.deg2rad(angle)
        self._detail_level = detail_level

    def plot(self, **kwargs):

        if np.atleast_1d(self._xo).shape[0] == 1:

            X, Y, Z = self._build_cone(
                self._xp, self._yp, self._zp, self._xo, self._yo, self._zo
            )

        else:

            X = []
            Y = []
            Z = []

            for xp, yp, zp, xo, yo, zo in zip(
                self._xp, self._yp, self._zp, self._xo, self._yo, self._zo
            ):

                x, y, z = self._build_cone(xp, yp, zp, xo, yo, zo)

                X.append(x)
                Y.append(y)
                Z.append(z)

        return ipv.plot_wireframe(X, Y, Z, **kwargs)

    def _build_cone(self, xp, yp, zp, xo, yo, zo):

        p0 = np.array([xo, yo, zo])

        vec = np.array([xp, yp, zp]) - p0

        height = np.linalg.norm(vec)

        vec = vec / height

        # make some vector not in the same direction as v

        not_v = np.array([1, 1, 0])

        if (vec == not_v).all():
            not_v = np.array([0, 1, 0])

        # make vector perpendicular to v

        n1 = np.cross(vec, not_v)
        # normalize n1

        n1 /= np.linalg.norm(n1)

        # make unit vector perpendicular to v and n1

        n2 = np.cross(vec, n1)

        # make unit vector perpendicular to v and n1

        n2 = np.cross(vec, n1)

        # surface ranges over t from 0 to length of axis and 0 to 2*pi

        t = np.linspace(0, height, self._detail_level)

        theta = np.linspace(0, 2 * np.pi, self._detail_level)

        # use meshgrid to make 2d arrays
        t, theta = np.meshgrid(t, theta)

        R = height * np.tan(self._angle)

        R = np.linspace(0, R, self._detail_level)

        # generate coordinates for surface
        X, Y, Z = [
            p0[i] + vec[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i]
            for i in [0, 1, 2]
        ]

        return X, Y, Z
