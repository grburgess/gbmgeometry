import ipyvolume as ipv
import numpy as np
import h5py

from gbmgeometry.utils.package_utils import get_path_of_data_file
from gbmgeometry.geometry import Sphere




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

        with h5py.File(get_path_of_data_file('countries.h5'), 'r') as f:

            xs = f['x'][()]
            ys = f['y'][()]
            zs = f['z'][()]
            
            ipv.pylab.plot(xs,ys,zs, color='black')


class Sol(Sphere):
    def __init__(self, x, y, z, ax=None, detail_level=50, color="#EF990F", **kwargs):
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
    def __init__(self, x, y, z, ax=None, detail_level=20, color="#68696A", **kwargs):
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
