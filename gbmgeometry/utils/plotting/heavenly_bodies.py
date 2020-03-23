import ipyvolume as ipv
import numpy as np

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


class Earth(Sphere):
    def __init__(self, ax=None, detail_level=100, color="#4E66DE", **kwargs):
        """
        The Planet Earth. Earth is at the origin

        :param ax: 
        :param detail_level: 
        :param color: 
        :returns: 
        :rtype: 

        """

        super(Earth, self).__init__(
            ax=ax,
            x=0,
            y=0,
            z=0,
            detail_level=detail_level,
            radius=6371.0,
            color=color,
            **kwargs
        )


class Sol(Sphere):
    def __init__(self, x, y, z, ax=None, detail_level=100, color="#EF990F", **kwargs):
        """
        The Sun. This is variable with respect to the satellite and Earth, so 
        coordinates have to be provided with at a single time or as an array

        :param ax: 
        :param detail_level: 
        :param color: 
        :returns: 
        :rtype: 

        """

        super(Sol, self).__init__(
            ax=ax,
            x=x,
            y=y,
            z=z,
            detail_level=detail_level,
            radius=696340.0,
            color=color,
            **kwargs
        )


class Moon(Sphere):
    def __init__(self, x, y, z, ax=None, detail_level=100, color="#68696A", **kwargs):
        """
        The Sun. This is variable with respect to the satellite and Earth, so 
        coordinates have to be provided with at a single time or as an array

        :param ax: 
        :param detail_level: 
        :param color: 
        :returns: 
        :rtype: 

        """

        super(Moon, self).__init__(
            ax=ax,
            x=x,
            y=y,
            z=z,
            detail_level=detail_level,
            radius=1731.1,
            color=color,
            **kwargs
        )


def _xyz(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def _sample_theta_phi(size, r):

    theta = np.arccos(1 - 2 * np.random.uniform(0.0, 1.0, size=size))
    phi = np.random.uniform(0, 2 * np.pi, size=size)

    x, y, z = _xyz(r, theta, phi)
    return x, y, z


class StarField(object):
    def __init__(self, n_stars=20, distance=1):

        self._x, self._y, self._z = _sample_theta_phi(n_stars, distance)

    def plot(self):

        return ipv.pylab.scatter(
            self._x, self._y, self._z, color="white", marker="sphere", size=.05, alpha=.2
        )
