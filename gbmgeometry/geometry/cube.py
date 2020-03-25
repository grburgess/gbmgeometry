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
    ):
        factor = 1.0

        if transform_matrix is not None:

            factor = 1
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

    def _create_Ploy3DCollection(
            self, *two_lines, i = None, side=None
    ):

        x, y, z = zip(two_lines[0], two_lines[1], two_lines[2], two_lines[3])
        if i is None:

            

            # print(x)

            X = np.array([__ for __ in iterable_to_chunks(x, 2)])
            Y = np.array([__ for __ in iterable_to_chunks(y, 2)])
            Z = np.array([__ for __ in iterable_to_chunks(z, 2)])

            return ipv.plot_surface(X, Y, Z, color=self._color)

        else:

            self._X[side,i] = np.array([__ for __ in iterable_to_chunks(x, 2)])
            self._Y[side,i] = np.array([__ for __ in iterable_to_chunks(y, 2)])
            self._Z[side,i] = np.array([__ for __ in iterable_to_chunks(z, 2)])

            if i >= len(self._x) -1:

                return ipv.plot_surface(self._X[side], self._Y[side], self._Z[side], color=self._color)
        
    @property
    def artists(self):
        return self._artists

    def plot(self):

        if np.atleast_1d(self._x).shape[0] == 1:

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

            self._X = np.zeros((6, len(self._x), 2, 2))
            self._Y = np.zeros((6, len(self._x), 2, 2))
            self._Z = np.zeros((6, len(self._x), 2, 2))

            for i, (x, y, z) in enumerate(zip(self._x, self._y, self._z)):

                (
                    front_top_right,
                    front_top_left,
                    front_bot_right,
                    front_bot_left,
                    rear_top_right,
                    rear_top_left,
                    rear_bot_right,
                    rear_bot_left,
                ) = self._assemble(x, y, z)

                top = self._create_Ploy3DCollection(
                     front_top_left, front_top_right, rear_top_left, rear_top_right,
                    i=i, side=0)

                bot = self._create_Ploy3DCollection(
                    front_bot_left, front_bot_right, rear_bot_left, rear_bot_right,
                    i=i, side=1)

                front = self._create_Ploy3DCollection(
                    front_bot_left, front_bot_right, front_top_left, front_top_right,
                    i=i, side=2)

                rear = self._create_Ploy3DCollection(
                    rear_top_left, rear_top_right, rear_bot_left, rear_bot_right,
                    i=i, side=3)

                left = self._create_Ploy3DCollection(
                    front_bot_left, front_top_left, rear_bot_left, rear_top_left,
                    i=i, side=4)

                right = self._create_Ploy3DCollection(
                    front_bot_right, front_top_right, rear_bot_right, rear_top_right,
                    i=i, side=5)

        self._artists.extend([top, bot, front, rear, left, right])        

    def _assemble(self, x, y, z):

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
            new_points = []
            for point in points:

                new_point = np.dot(self._transform_matrix, point)
                new_point = self._translate(new_point)
                #                new_point = self._translate(self._quaternion.rotatePoint(point))

                # print(f"before {point} after {new_point}")

                new_points.append(new_point)

            # (
            #     front_top_right,
            #     front_top_left,
            #     front_bot_right,
            #     front_bot_left,
            #     rear_top_right,
            #     rear_top_left,
            #     rear_bot_right,
            #     rear_bot_left,
            # ) = new_points

            return new_points

        else:

            return points

    def _translate(self, point):

        return np.array(
            [
                point[0] + self._sc_pos[0],
                point[1] + self._sc_pos[1],
                point[2] + self._sc_pos[2],
            ]
        )


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
