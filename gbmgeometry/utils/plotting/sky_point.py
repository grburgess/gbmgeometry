import ipyvolume as ipv
from astropy.coordinates import SkyCoord
import numpy as np


class SkyPoint(object):
    def __init__(self, ra, dec, color="#33CF64"):
        self._position = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
        self._color = color
        self._x, self._y, self._x = self._position.cartesian.xyz.to("km").value

    def plot(self, sx, sy, sz):

        ipv.pylab.plot(
            np.array([sx, sx + self._x]),
            np.array([sy, sy + self._y]),
            np.array([sz, sz + self._z]),
            color=self._color,
        )
