import numpy as np
from astropy.coordinates import frame_transform_graph
from astropy.coordinates import BaseCoordinateFrame, FrameAttribute, TimeFrameAttribute, RepresentationMapping
import astropy.coordinates as coord
import astropy.units as u


class GBMFrame(BaseCoordinateFrame):
    """
    
    Fermi GBM Frame

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
  
    """
    default_representation = coord.SphericalRepresentation

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping(reprname='lon', framename='Az', defaultunit=u.degree),
                      RepresentationMapping(reprname='lat', framename='Zen', defaultunit=u.degree),
                      RepresentationMapping(reprname='distance', framename='DIST', defaultunit=None)],
        'unitspherical': [RepresentationMapping(reprname='lon', framename='Az', defaultunit=u.degree),
                          RepresentationMapping(reprname='lat', framename='Zen', defaultunit=u.degree)],
        'cartesian': [RepresentationMapping(reprname='x', framename='SCX'),
                      RepresentationMapping(reprname='y', framename='SCY'),
                      RepresentationMapping(reprname='z', framename='SCZ')]
    }

    # Specify frame attributes required to fully specify the frame
    # location = FrameAttribute(default=None)
    quaternion = FrameAttribute(default=None)

    # equinox = TimeFrameAttribute(default='J2000')

    @staticmethod
    def _set_quaternion(quaternion):
        sc_matrix = np.zeros((3, 3))

        sc_matrix[0, 0] = (quaternion[0] ** 2 - quaternion[1] ** 2 - quaternion[2] ** 2 + quaternion[3] ** 2)
        sc_matrix[0, 1] = 2.0 * (quaternion[0] * quaternion[1] + quaternion[3] * quaternion[2])
        sc_matrix[0, 2] = 2.0 * (quaternion[0] * quaternion[2] - quaternion[3] * quaternion[1])
        sc_matrix[1, 0] = 2.0 * (quaternion[0] * quaternion[1] - quaternion[3] * quaternion[2])
        sc_matrix[1, 1] = (-quaternion[0] ** 2 + quaternion[1] ** 2 - quaternion[2] ** 2 + quaternion[3] ** 2)
        sc_matrix[1, 2] = 2.0 * (quaternion[1] * quaternion[2] + quaternion[3] * quaternion[0])
        sc_matrix[2, 0] = 2.0 * (quaternion[0] * quaternion[2] + quaternion[3] * quaternion[1])
        sc_matrix[2, 1] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[3] * quaternion[0])
        sc_matrix[2, 2] = (-quaternion[0] ** 2 - quaternion[1] ** 2 + quaternion[2] ** 2 + quaternion[3] ** 2)

        return sc_matrix


@frame_transform_graph.transform(coord.FunctionTransform, GBMFrame, coord.ICRS)
def gbm_to_j2000(gbm_coord, j2000_frame):
    """ Compute the transformation from heliocentric Sgr coordinates to
        spherical Galactic.
    """

    sc_matrix = gbm_coord._set_quaternion(gbm_coord.quaternion)

    # X,Y,Z = gbm_coord.cartesian

    pos = gbm_coord.cartesian.xyz.value

    X0 = np.dot(sc_matrix[:, 0], pos)
    X1 = np.dot(sc_matrix[:, 1], pos)
    X2 = np.clip(np.dot(sc_matrix[:, 2], pos), -1., 1.)

    dec = np.arcsin(X2)

    idx = np.logical_and(np.abs(X0) < 1E-6, np.abs(X1) < 1E-6)

    ra = np.zeros_like(dec)

    ra[~idx] = np.arctan2(X1, X0) % (2 * np.pi)

    return coord.ICRS(ra=ra * u.radian, dec=dec * u.radian)


@frame_transform_graph.transform(coord.FunctionTransform, coord.ICRS, GBMFrame)
def j2000_to_gbm(j2000_frame, gbm_coord):
    """ Compute the transformation from heliocentric Sgr coordinates to
        spherical Galactic.
    """

    sc_matrix = gbm_coord._set_quaternion(gbm_coord.quaternion)

    pos = j2000_frame.cartesian.xyz.value

    X0 = np.dot(sc_matrix[0, :], pos)
    X1 = np.dot(sc_matrix[1, :], pos)
    X2 = np.dot(sc_matrix[2, :], pos)

    el = np.pi / 2. - np.arccos(X2)  # convert to proper frame

    idx = np.logical_and(np.abs(X0) < 1E-6, np.abs(X1) < 1E-6)

    az = np.zeros_like(el)

    az[~idx] = np.arctan2(X1, X0) % (2 * np.pi)

    return GBMFrame(Az=az * u.radian, Zen=el * u.radian, quaternion=gbm_coord.quaternion)
