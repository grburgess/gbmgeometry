import numpy as np
import collections



class Volume(object):
    def __init__(self, x_origin, y_origin, z_origin, height, x_width, y_width, color='grey'):
        self._center = (x_origin, y_origin, z_origin)

        self._color = color


        self._build_cube(origin=(x_origin - x_width / 2., z_origin - height / 2., y_origin - y_width / 2.),
                         width=x_width,
                         depth=y_width,
                         height=height)




    def _build_cube(self,origin=None, width=1,height=1, depth=1):



        self._planes = collections.OrderedDict()
        for plane in ['+x', '-x', '+y', '-y', '+z', '-z']:

            self._planes[plane] = None



        u, v, w = (0, 0, 0) if origin is None else origin


        grids = []

        for plane in self._planes.keys():


            if '-z' in plane:

                this_grid = self._grid('xy', (u, w), width, depth, v)

                self._planes[plane] = this_grid

                grids.extend(this_grid)

            if '+z' in plane:

                this_grid = self._grid('xy', (u, w), width, depth, v + height)

                self._planes[plane] = this_grid

                grids.extend(this_grid)


            if '-y' in plane:


                this_grid = self._grid('xz', (u, v), width, height, w)

                self._planes[plane] = this_grid

                grids.extend(this_grid)

            if '+y' in plane:

                this_grid = self._grid('xz', (u, v), width, height, w + depth)

                self._planes[plane] = this_grid

                grids.extend(this_grid)

            if '-x' in plane:

                this_grid = self._grid('yz', (w, v), depth, height, u)

                self._planes[plane] = this_grid

                grids.extend(this_grid)

            if '+x' in plane:

                this_grid = self._grid('yz', (w, v), depth, height, u + width)

                self._planes[plane] = this_grid

                grids.extend(this_grid)

        self._quads = np.array(grids)


    @staticmethod
    def _grid(plane='xy',origin=None, width=1 ,height=1, depth=0):

        u, v = (0, 0) if origin is None else origin

        # hard code this
        width_segments = 1
        height_segments = 1

        w_x, h_y = width, height

        quads = []

        for i in range(width_segments):
            for j in range(height_segments):
                quads.append(
                    Volume._quad(plane, (i * w_x + u, j * h_y + v), w_x, h_y, depth))

        return np.array(quads)

    @staticmethod
    def _quad(plane='xy', origin=None, width=1, height=1, depth=0):
        u, v = (0, 0) if origin is None else origin

        plane = plane.lower()
        if plane == 'xy':
            vertices = ((u, v, depth),
                        (u + width, v, depth),
                        (u + width, v + height, depth),
                        (u, v + height, depth))
        elif plane == 'xz':
            vertices = ((u, depth, v),
                        (u + width, depth, v),
                        (u + width, depth, v + height),
                        (u, depth, v + height))
        elif plane == 'yz':
            vertices = ((depth, u, v),
                        (depth, u + width, v),
                        (depth, u + width, v + height),
                        (depth, u, v + height))
        else:
            raise ValueError('"{0}" is not a supported plane!'.format(plane))

        return np.array(vertices)



    @property
    def quads(self):
        return self._quads

    @property
    def center(self):
        return self._center

    def plot(self, ax):
        collection = Poly3DCollection(self._quads,
                                      facecolors=self._color,
                                      alpha=.1)

        collection.set_facecolor((0, 0, 1, .1))

        ax.add_collection3d(collection)

    def intersect(self, ray):
        pass

