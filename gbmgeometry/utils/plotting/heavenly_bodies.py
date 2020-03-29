import ipyvolume as ipv
import numpy as np
import h5py
import PIL.Image as pil_image
from gbmgeometry.utils.package_utils import get_path_of_data_file
from gbmgeometry.geometry import Sphere
from gbmgeometry.geometry.cirs_to_gcrs import cirs_to_gcrs

_earth_img = dict(
    day=get_path_of_data_file("earth_day.jpg"),
    night=get_path_of_data_file("earth_night.jpg"),
    midnight=get_path_of_data_file("earth_midnight.jpg"),
    thats_no_moon=get_path_of_data_file("thats_no_moon.jpg"),
    i_have_a_bad_feeling_about_this=get_path_of_data_file("thats_no_moon_big.jpg"),
)


class Earth(Sphere):
    def __init__(
        self,
        ax=None,
        detail_level=50,
        color="#040C6A",
        earth_time="day",
        realistic=True,
        astro_time=None,
        **kwargs,
    ):
        """
        The Planet Earth. Earth is at the origin

        :param ax: 
        :param detail_level: 
        :param color: 
        :returns: 
        :rtype: 

        """

        transform_matrix = None

        if realistic:

            assert astro_time is not None

            assert (
                earth_time in _earth_img
            ), "oops, please choose, day, night, or midnight"

            image = pil_image.open(_earth_img[earth_time])

            # compute the transformtion matrix

            if np.atleast_1d(astro_time).shape[0] > 1:

                transform_matrix = np.zeros((len(astro_time), 3, 3))

                for i in range(len(astro_time)):
                    transform_matrix[i] = cirs_to_gcrs(astro_time[i])

            else:

                transform_matrix = cirs_to_gcrs(astro_time)

        else:

            image = None

        super(Earth, self).__init__(
            ax=ax,
            x=0,
            y=0,
            z=0,
            detail_level=detail_level,
            radius=6371.0,
            color=color,
            image=image,
            transform_matrix=transform_matrix,
            **kwargs,
        )

        if not realistic:

            with h5py.File(get_path_of_data_file("countries.h5"), "r") as f:

                xs = f["x"][()]
                ys = f["y"][()]
                zs = f["z"][()]

                ipv.pylab.plot(xs, ys, zs, color="black")


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
            **kwargs,
        )


class Moon(Sphere):
    def __init__(
        self,
        x,
        y,
        z,
        ax=None,
        detail_level=30,
        color="#68696A",
        realistic=False,
        **kwargs,
    ):
        """
        The Sun. This is variable with respect to the satellite and Earth, so 
        coordinates have to be provided with at a single time or as an array

        :param ax: 
        :param detail_level: 
        :param color: 
        :returns: 
        :rtype: 

        """

        # it is crazy slow to animate and image!

        if realistic:

            image = pil_image.open(_earth_img["thats_no_moon"])

        else:
            image = None

        super(Moon, self).__init__(
            ax=ax,
            x=x,
            y=y,
            z=z,
            detail_level=detail_level,
            radius=1731.1,
            color=color,
            image=image,
            **kwargs,
        )


def _xyz(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def _sample_theta_phi(size, r):

    theta = np.arccos(1 - 2 * np.random.uniform(0.0, 1.0, size=size))
    phi = np.random.uniform(0, 2 * np.pi, size=size)
    r = np.random.uniform(r * 0.95, r * 1.1, size=size)

    x, y, z = _xyz(r, theta, phi)
    return x, y, z


class StarField(object):
    def __init__(self, n_stars=20, distance=1):
        """
        a star field

        :param n_stars: 
        :param distance: 
        :returns: 
        :rtype: 

        """

        self._x, self._y, self._z = _sample_theta_phi(n_stars, distance)

    def plot(self):

        return ipv.pylab.scatter(
            self._x,
            self._y,
            self._z,
            color="white",
            marker="sphere",
            size=0.07,
            alpha=0.2,
        )
