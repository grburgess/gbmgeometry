import ipyvolume as ipv
from astropy.coordinates import SkyCoord
import numpy as np
import healpy as hp
from gbmgeometry.utils.array_to_cmap import array_to_cmap


class SkyPoint(object):
    def __init__(
        self, ra, dec, color="#33CF64", alpha=1.0, distance=1e10, as_point=False
    ):
        """

        A point on the sky

        :param ra: 
        :param dec: 
        :param color: 
        :param distance: 
        :returns: 
        :rtype: 

        """

        self._position = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
        self._color = color

        if as_point:
            self._x, self._y, self._z = self._position.cartesian.xyz.value

        else:

            self._x, self._y, self._z = self._position.cartesian.xyz.value * distance

        self._as_point = as_point
        self._distance = distance

    def plot(self, sx, sy, sz, distance=None):

        if self._as_point:

            return ipv.pylab.scatter(
                np.array([sx + self._x * distance]),
                np.array([sy + self._y * distance]),
                np.array([sz + self._z * distance]),
                color=self._color,
                marker="circle_2d",
                size=0.5,
            )

        else:

            return ipv.pylab.plot(
                np.array([sx, sx + self._x]),
                np.array([sy, sy + self._y]),
                np.array([sz, sz + self._z]),
                color=self._color,
            )


def balrog_to_skypoints(
    healpix_map_file, cmap="viridis", new_nside=None, as_point=False
):

    # collect the nside of the map

    healpix_map = hp.read_map(healpix_map_file)

    if new_nside is not None:

        healpix_map = hp.ud_grade(healpix_map, new_nside)

    nside = hp.get_nside(healpix_map)

    # get the colors of the rays based of their values

    _, colors = array_to_cmap(healpix_map, cmap=cmap, use_log=False)

    # now go thru all the points on the sphere

    sky_points = []

    for idx in healpix_map.nonzero()[0]:

        ra, dec = pix_to_sky(idx, nside)

        # mark the color

        color = colors[idx]

        sky_points.append(SkyPoint(ra, dec, color=color, alpha=0.7, as_point=as_point))

    return sky_points


def pix_to_sky(idx, nside):
    """Convert the pixels corresponding to the input indexes to sky coordinates (RA, Dec)"""

    theta, phi = hp.pix2ang(nside, idx)

    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    return ra, dec
