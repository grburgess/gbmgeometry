from itertools import chain, repeat
import numpy as np
import ipyvolume as ipv
from gbmgeometry.geometry import Quaternion


class Cube(object):
    def __init__(
        self,
        x=0,
        y=0,
        z=0,
        x_width=1.0,
        y_width=1.0,
        z_width=1.0,
        color="red",
        transform_matrix=None,
        sc_pos=None,
        quarternion=None,
        scaling_factor=1.0,
    ):
        """
        A 3D cube that can be animated

        :param x: 
        :param y: 
        :param z: 
        :param x_width: 
        :param y_width: 
        :param z_width: 
        :param color: 
        :param transform_matrix: 
        :param sc_pos: 
        :param quarternion: 
        :param scaling_factor: 
        :returns: 
        :rtype: 

        """

        factor = 1.0

        if transform_matrix is not None:

            factor = scaling_factor
            x *= factor
            y *= factor
            z *= factor
            # also convert to

        self._x = x - 0.5 * x_width
        self._y = y - 0.5 * y_width
        self._z = z - 0.5 * z_width
        self._sc_pos = sc_pos

        # print(f"cube {x} {y} {z}")
        # print(f"corner {self._x} {self._y} {self._z}")

        self._height = y_width
        self._width = x_width
        self._depth = z_width

        self._color = color
        self._artists = []
        self._transform_matrix = transform_matrix

        self._is_vector = False
        self._is_vector_origin = False

        # check if we have vector of positions

        if (self._transform_matrix is not None) and len(
            self._transform_matrix.shape
        ) > 2:

            # this is the case if the transform is changing with time
            # or the sc_pos is changing with time

            assert (
                self._sc_pos.shape[0] == self._transform_matrix.shape[0]
            ), "the sc_pos and trans matrix are not the same lenght"

            self._is_vector = True
            self._vector_length = self._transform_matrix.shape[0]

        elif np.atleast_1d(self._x).shape[0] > 1:

            # this is the case that the origin is changing with time.

            self._is_vector = True
            self._is_vector_origin = True
            self._vector_length = self._x.shape[0]

    def _create_Ploy3DCollection(self, line1, line2, line3, line4, i=None, side=None):

        x, y, z = zip(line1, line2, line3, line4)
        if not self._is_vector:

            # print(x)

            X = np.array([__ for __ in iterable_to_chunks(x, 2)])
            Y = np.array([__ for __ in iterable_to_chunks(y, 2)])
            Z = np.array([__ for __ in iterable_to_chunks(z, 2)])

            return ipv.plot_surface(X, Y, Z, color=self._color)

        else:

            self._X[side, i] = np.array([__ for __ in iterable_to_chunks(x, 2)])
            self._Y[side, i] = np.array([__ for __ in iterable_to_chunks(y, 2)])
            self._Z[side, i] = np.array([__ for __ in iterable_to_chunks(z, 2)])

            # on the last iteration we make the plot and return it

            if i >= self._vector_length - 1:

                return ipv.plot_surface(
                    self._X[side], self._Y[side], self._Z[side], color=self._color
                )

    @property
    def artists(self):
        return self._artists

    def plot(self):

        if not self._is_vector:

            (
                front_top_right,
                front_top_left,
                front_bot_right,
                front_bot_left,
                rear_top_right,
                rear_top_left,
                rear_bot_right,
                rear_bot_left,
            ) = self._assemble(self._x, self._y, self._z)

            top = self._create_Ploy3DCollection(
                front_top_left, front_top_right, rear_top_left, rear_top_right,
            )

            bot = self._create_Ploy3DCollection(
                front_bot_left, front_bot_right, rear_bot_left, rear_bot_right,
            )

            front = self._create_Ploy3DCollection(
                front_bot_left, front_bot_right, front_top_left, front_top_right,
            )

            rear = self._create_Ploy3DCollection(
                rear_top_left, rear_top_right, rear_bot_left, rear_bot_right,
            )

            left = self._create_Ploy3DCollection(
                front_bot_left, front_top_left, rear_bot_left, rear_top_left,
            )

            right = self._create_Ploy3DCollection(
                front_bot_right, front_top_right, rear_bot_right, rear_top_right,
            )

        else:

            self._X = np.zeros((6, self._vector_length, 2, 2))
            self._Y = np.zeros((6, self._vector_length, 2, 2))
            self._Z = np.zeros((6, self._vector_length, 2, 2))

            for i in range(self._vector_length):

                if self._is_vector_origin:
                    x = self._x[i]
                    y = self._y[i]
                    z = self._z[i]

                else:

                    x, y, z = self._x, self._y, self._z

                (
                    front_top_right,
                    front_top_left,
                    front_bot_right,
                    front_bot_left,
                    rear_top_right,
                    rear_top_left,
                    rear_bot_right,
                    rear_bot_left,
                ) = self._assemble(x, y, z, i)

                top = self._create_Ploy3DCollection(
                    front_top_left,
                    front_top_right,
                    rear_top_left,
                    rear_top_right,
                    i=i,
                    side=0,
                )

                bot = self._create_Ploy3DCollection(
                    front_bot_left,
                    front_bot_right,
                    rear_bot_left,
                    rear_bot_right,
                    i=i,
                    side=1,
                )

                front = self._create_Ploy3DCollection(
                    front_bot_left,
                    front_bot_right,
                    front_top_left,
                    front_top_right,
                    i=i,
                    side=2,
                )

                rear = self._create_Ploy3DCollection(
                    rear_top_left,
                    rear_top_right,
                    rear_bot_left,
                    rear_bot_right,
                    i=i,
                    side=3,
                )

                left = self._create_Ploy3DCollection(
                    front_bot_left,
                    front_top_left,
                    rear_bot_left,
                    rear_top_left,
                    i=i,
                    side=4,
                )

                right = self._create_Ploy3DCollection(
                    front_bot_right,
                    front_top_right,
                    rear_bot_right,
                    rear_top_right,
                    i=i,
                    side=5,
                )

        self._artists.extend([top, bot, front, rear, left, right])

    def _assemble(self, x, y, z, i=None):

        front_bot_left = np.array([x, y, z])
        front_bot_right = np.array([x + self._width, y, z])
        front_top_left = np.array([x, y + self._height, z])
        front_top_right = np.array([x + self._width, y + self._height, z])

        rear_bot_left = np.array([x, y, z + self._depth])
        rear_bot_right = np.array([x + self._width, y, z + self._depth])
        rear_top_left = np.array([x, y + self._height, z + self._depth])
        rear_top_right = np.array([x + self._width, y + self._height, z + self._depth,])

        points = [
            front_top_right,
            front_top_left,
            front_bot_right,
            front_bot_left,
            rear_top_right,
            rear_top_left,
            rear_bot_right,
            rear_bot_left,
        ]

        if self._transform_matrix is not None:

            # we now rotate the point FIRST
            # and then translate the point

            new_points = []

            for point in points:

                new_point = self._rotate(point, i)

                new_point = self._translate(new_point, i)

                # print(f"before {point} after {new_point}")

                new_points.append(new_point)

            return new_points

        else:

            return points

    def _translate(self, point, i=None):

        if i is not None:

            return point + self._sc_pos[i]

        else:

            return point + self._sc_pos

    def _rotate(self, point, i=None):

        if i is not None:

            return np.dot(self._transform_matrix[i].T, point)

        else:

            return np.dot(self._transform_matrix.T, point)


def iterable_to_chunks(iterable, size, fill=None):
    """
    Split a list to chunks
    iterable : an iterable object, e.g. generator, list, array, etc.
    size : chunk size, positive integer
    fill : padding values.
    Example :
      tochunks('abcdefg', 3, 'x')
    Output :
      ('a','b','c'), ('d','e','f'), ('g','x','x')
    reference :
      http://stackoverflow.com/a/312644
    """
    return zip(*[chain(iterable, repeat(fill, size - 1))] * size)
