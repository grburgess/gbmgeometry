import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from astropy.coordinates import BaseCoordinateFrame, Attribute, RepresentationMapping
from astropy.coordinates import frame_transform_graph


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
        "spherical": [
            RepresentationMapping(
                reprname="lon", framename="lon", defaultunit=u.degree
            ),
            RepresentationMapping(
                reprname="lat", framename="lat", defaultunit=u.degree
            ),
            RepresentationMapping(
                reprname="distance", framename="DIST", defaultunit=None
            ),
        ],
        "unitspherical": [
            RepresentationMapping(
                reprname="lon", framename="lon", defaultunit=u.degree
            ),
            RepresentationMapping(
                reprname="lat", framename="lat", defaultunit=u.degree
            ),
        ],
        "cartesian": [
            RepresentationMapping(reprname="x", framename="SCX"),
            RepresentationMapping(reprname="y", framename="SCY"),
            RepresentationMapping(reprname="z", framename="SCZ"),
        ],
    }

    # Specify frame attributes required to fully specify the frame
    sc_pos_X = Attribute(default=None)
    sc_pos_Y = Attribute(default=None)
    sc_pos_Z = Attribute(default=None)

    quaternion_1 = Attribute(default=None)
    quaternion_2 = Attribute(default=None)
    quaternion_3 = Attribute(default=None)
    quaternion_4 = Attribute(default=None)

    # equinox = TimeFrameAttribute(default='J2000')

    @staticmethod
    def _set_quaternion(q1, q2, q3, q4):
        sc_matrix = np.zeros((3, 3))

        sc_matrix[0, 0] = q1 ** 2 - q2 ** 2 - q3 ** 2 + q4 ** 2
        sc_matrix[0, 1] = 2.0 * (q1 * q2 + q4 * q3)
        sc_matrix[0, 2] = 2.0 * (q1 * q3 - q4 * q2)
        sc_matrix[1, 0] = 2.0 * (q1 * q2 - q4 * q3)
        sc_matrix[1, 1] = -(q1 ** 2) + q2 ** 2 - q3 ** 2 + q4 ** 2
        sc_matrix[1, 2] = 2.0 * (q2 * q3 + q4 * q1)
        sc_matrix[2, 0] = 2.0 * (q1 * q3 + q4 * q2)
        sc_matrix[2, 1] = 2.0 * (q2 * q3 - q4 * q1)
        sc_matrix[2, 2] = -(q1 ** 2) - q2 ** 2 + q3 ** 2 + q4 ** 2

        return sc_matrix


@frame_transform_graph.transform(coord.FunctionTransform, GBMFrame, coord.ICRS)
def gbm_to_j2000(gbm_coord, j2000_frame):
    """ Compute the transformation from heliocentric Sgr coordinates to
        spherical Galactic.
    """

    sc_matrix = gbm_coord._set_quaternion(
        gbm_coord.quaternion_1,
        gbm_coord.quaternion_2,
        gbm_coord.quaternion_3,
        gbm_coord.quaternion_4,
    )

    # X,Y,Z = gbm_coord.cartesian

    pos = gbm_coord.cartesian.xyz.value

    X0 = np.dot(sc_matrix[:, 0], pos)
    X1 = np.dot(sc_matrix[:, 1], pos)
    X2 = np.clip(np.dot(sc_matrix[:, 2], pos), -1.0, 1.0)

    # dec = np.arcsin(X2)

    dec = np.pi / 2.0 - np.arccos(X2)

    idx = np.logical_and(np.abs(X0) < 1e-6, np.abs(X1) < 1e-6)

    ra = np.zeros_like(dec)

    ra[~idx] = np.arctan2(X1, X0) % (2 * np.pi)

    return coord.ICRS(ra=ra * u.radian, dec=dec * u.radian)


@frame_transform_graph.transform(coord.FunctionTransform, coord.ICRS, GBMFrame)
def j2000_to_gbm(j2000_frame, gbm_coord):
    """ Compute the transformation from heliocentric Sgr coordinates to
        spherical Galactic.
    """

    sc_matrix = gbm_coord._set_quaternion(
        gbm_coord.quaternion_1,
        gbm_coord.quaternion_2,
        gbm_coord.quaternion_3,
        gbm_coord.quaternion_4,
    )

    pos = j2000_frame.cartesian.xyz.value

    X0 = np.dot(sc_matrix[0, :], pos)
    X1 = np.dot(sc_matrix[1, :], pos)
    X2 = np.clip(np.dot(sc_matrix[2, :], pos), -1.0, 1.0)
    el = np.pi / 2.0 - np.arccos(X2)  # convert to proper frame

    idx = np.logical_and(np.abs(X0) < 1e-6, np.abs(X1) < 1e-6)

    az = np.zeros_like(el)

    az[~idx] = np.arctan2(X1, X0) % (2 * np.pi)

    az[np.rad2deg(el) == 90.0] = 0.0

    return GBMFrame(
        lon=az * u.radian,
        lat=el * u.radian,
        quaternion_1=gbm_coord.quaternion_1,
        quaternion_2=gbm_coord.quaternion_2,
        quaternion_3=gbm_coord.quaternion_3,
        quaternion_4=gbm_coord.quaternion_4,
    )
