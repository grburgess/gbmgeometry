### Functions to add small cylinders as detector crystals ###

import numpy as np
import math


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    :param axis: Axis vector to rotate around
    :param theta: Angle to rotate
    :return: Rotation matrix
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def add_rotated_cylinder(
    ax,
    height=30,
    radius=15,
    x_center=0,
    y_center=0,
    z_center=0,
    color="magenta",
    alpha=0.5,
    theta=0,
    phi=0,
    ceiling=True,
    floor=True,
    label=None,
):
    """
    Returns the unit cylinder that corresponds to the curve r.
    :param ax: ax of figure to add the cylinder
    :param height: Height of cylinder
    :param radius: Radius of cylinder
    :param x_center: x-coord of center of cylinder
    :param y_center: y-coord of center of cylinder
    :param z_center: z-coord of center of cylinder
    :param color: color of cylinder
    :param alpha: alpha of cylinder
    :param theta: theta angle of cylinder orientation
    :param phi: phi angle of cylinder orientation
    :param ceiling: add a ceiling to cylinder?
    :param floor: add a floor to cylinder?
    :param label: Label for det
    :return:
    """

    r = np.ones(10)
    n = 20
    # ensure that r is a column vector
    r = np.atleast_2d(r)
    r_rows, r_cols = r.shape

    if r_cols > r_rows:
        r = r.T

    # find points along x and y axes
    points = np.linspace(0, 2 * np.pi, n + 1)
    x = np.cos(points) * r * radius  # +x_center
    y = np.sin(points) * r * radius  # +y_center

    # find points along z axis
    rpoints = np.atleast_2d(np.linspace(0, 1, len(r)))
    z = np.ones((1, n + 1)) * rpoints.T * height  # +z_center

    axis1 = [0, 1.0, 0]
    axis2 = [0, 0, 1.0]
    xr, yr, zr = np.dot(
        rotation_matrix(axis2, phi),
        np.transpose(
            np.dot(
                rotation_matrix(axis1, theta),
                np.transpose(np.array([x, y, z]), (1, 0, 2)),
            ),
            (1, 0, 2),
        ),
    )
    if label is not None:
        surf = ax.plot_surface(
            xr + x_center,
            yr + y_center,
            zr + z_center,
            color=color,
            alpha=alpha,
            label=label,
        )
    else:
        surf = ax.plot_surface(
            xr + x_center, yr + y_center, zr + z_center, color=color, alpha=alpha
        )

    # surf._facecolors2d = surf._facecolors3d
    # surf._edgecolors2d = surf._edgecolors3d

    # Ceiling and floor
    if floor:
        R = np.linspace(0, radius, 10)
        h = 0  # z_center
        u = np.linspace(0, 2 * np.pi, 100)

        x = np.outer(R, np.cos(u))  # +x_center
        y = np.outer(R, np.sin(u))  # +y_center

        xr, yr, hr = np.dot(
            rotation_matrix(axis2, phi),
            np.transpose(
                np.dot(
                    rotation_matrix(axis1, theta),
                    np.transpose(np.array([x, y, h * np.ones_like(x)]), (1, 0, 2)),
                ),
                (1, 0, 2),
            ),
        )

        surf = ax.plot_surface(
            xr + x_center, yr + y_center, hr + z_center, color=color, alpha=alpha
        )
        # surf._facecolors2d = surf._facecolors3d
        # surf._edgecolors2d = surf._edgecolors3d

    # Celling and floor
    if ceiling:
        R = np.linspace(0, radius, 10)
        h = height  # z_center+height
        u = np.linspace(0, 2 * np.pi, 100)

        x = np.outer(R, np.cos(u))  # +x_center
        y = np.outer(R, np.sin(u))  # +y_center

        xr, yr, hr = np.dot(
            rotation_matrix(axis2, phi),
            np.transpose(
                np.dot(
                    rotation_matrix(axis1, theta),
                    np.transpose(np.array([x, y, h * np.ones_like(x)]), (1, 0, 2)),
                ),
                (1, 0, 2),
            ),
        )

        surf = ax.plot_surface(
            xr + x_center, yr + y_center, hr + z_center, color=color, alpha=alpha
        )
        # surf._facecolors2d = surf._facecolors3d
        # surf._edgecolors2d = surf._edgecolors3d
