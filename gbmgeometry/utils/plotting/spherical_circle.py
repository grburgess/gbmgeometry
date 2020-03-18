import numpy as np
from astropy.visualization.wcsaxes.patches import _rotate_polygon
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import astropy.units as u


class SphericalCircle(PathPatch):
    # created from the astropy.visualization.wcsaxes.patches.SphericalCircle class
    # changed to path from polygon to create disjointed parts
    """
    Create a patch representing a spherical circle - that is, a circle that is
    formed of all the points that are within a certain angle of the central
    coordinates on a sphere. Here we assume that latitude goes from -90 to +90

    This class is needed in cases where the user wants to add a circular patch
    to a celestial image, since otherwise the circle will be distorted, because
    a fixed interval in longitude corresponds to a different angle on the sky
    depending on the latitude.

    Parameters
    ----------
    center : tuple or `~astropy.units.Quantity`
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements.
    radius : `~astropy.units.Quantity`
        The radius of the circle
    resolution : int, optional
        The number of points that make up the circle - increase this to get a
        smoother circle.
    vertex_unit : `~astropy.units.Unit`
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """

    def __init__(self, ra, dec, radius, resolution=100, vertex_unit=u.degree, **kwargs):

        # Extract longitude/latitude, either from a tuple of two quantities, or
        # a single 2-element Quantity.

        longitude, latitude = ra, dec

        # #longitude values restricted on domain of (-180,180]
        # if longitude.to_value(u.deg) > 180. :
        # 	longitude = -360. * u.deg + longitude.to(u.deg)

        # Start off by generating the circle around the North pole
        lon = np.linspace(0.0, 2 * np.pi, resolution + 1)[:-1] * u.radian
        lat = np.repeat(0.5 * np.pi - np.deg2rad(radius), resolution) * u.radian

        lon, lat = _rotate_polygon(lon, lat, longitude, latitude)

        # Extract new longitude/latitude in the requested units
        lon = lon.to_value(vertex_unit)
        lat = lat.to_value(vertex_unit)
        # Create polygon vertices
        vertices = np.array([lon, lat]).transpose()

        # split path into two sections if circle crosses -180, 180 bounds
        codes = []
        last = (4000.4 * u.degree).to_value(
            vertex_unit
        )  # 400.4 is a random number large enough so first element is "MOVETO"
        for v in vertices:
            if np.absolute(v[0] - last) > (300 * u.degree).to_value(vertex_unit):
                codes.append(Path.MOVETO)
            else:
                codes.append(Path.LINETO)
            last = v[0]

        circle_path = Path(vertices, codes)

        super().__init__(circle_path, **kwargs)
