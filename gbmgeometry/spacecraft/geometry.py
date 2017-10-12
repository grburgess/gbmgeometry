import numpy as np
import collections
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from sympy import Plane, Point3D, Line3D





class Ray(object):

    R = 400. # cm


    def __init__(self,detector, point_source, probability=None):


        self._detector = detector

        self._probability = probability

        self._point_source = point_source

        self._calculate_ray_origin()

        self._sympy_line = Line3D( Point3D(self._origin), Point3D(self._detector.mount_point))


    def _calculate_ray_origin(self):


        theta = np.deg2rad(self._point_source.Zen.value)
        phi = np.deg2rad(self._point_source.Az.value)

        x = Ray.R * np.cos(phi) * np.sin(theta)
        y = Ray.R * np.sin(phi) * np.sin(theta)
        z = Ray.R * np.cos(theta)


        self._origin = np.array([x,y,z])


    def plot(self,ax):


        if self._probability is None:

            ax.plot([self.ray_origin[0],self.detector_origin[0]],
                    [self.ray_origin[1], self.detector_origin[1]],
                    [self.ray_origin[2], self.detector_origin[2]]


                    )


        else:

            raise NotImplementedError()

    @property
    def detector_name(self):

        return self._detector.name

    @property
    def probability(self):

        return self._probability

    @property
    def detector_origin(self):

        return self._detector.mount_point



    @property
    def ray_origin(self):

        return self._origin

    @property
    def sympy_line(self):

        return self._sympy_line



class Surface(object):

    def __init__(self, name, vertices):

        self._vertices = vertices

        self._name = name

        if '+' in name:

            self._sign = 1.

        elif '-' in name:

            self._sign = -1.

        else:

            raise RuntimeError('the plane name is wrong')


        self._calculate_origin()

    def _calculate_origin(self):
        """
        
        compute the origin of the plane
        
        
        :return: 
        """


        if 'x' in self._name:


            assert len(np.unique(self._vertices[:,0])) == 1, 'vertices are wrong!'

            self._normal = np.array([1,0,0]) * self._sign

            x_origin = self._vertices[0,0]

            y_origin = (min(self._vertices[:,1]) + max(self._vertices[:,1]))/2.

            z_origin = (min(self._vertices[:, 2]) + max(self._vertices[:, 2])) / 2.

            #self._edges = [self.]


        elif 'y' in self._name:

            assert len(np.unique(self._vertices[:, 1])) == 1, 'vertices are wrong!'

            self._normal = np.array([0, 1., 0]) * self._sign

            x_origin = (min(self._vertices[:, 0]) + max(self._vertices[:, 0])) / 2.

            y_origin = self._vertices[0, 1]

            z_origin = (min(self._vertices[:, 2]) + max(self._vertices[:, 2])) / 2.


        elif 'z' in self._name:

            assert len(np.unique(self._vertices[:, 2])) == 1, 'vertices are wrong!'

            self._normal = np.array([0, 0, 1]) * self._sign

            x_origin = (min(self._vertices[:, 0]) + max(self._vertices[:, 0])) / 2.

            y_origin = (min(self._vertices[:, 1]) + max(self._vertices[:, 1])) / 2.

            z_origin = self._vertices[0, 2]


        self._xmax = self._vertices[:,0].max()
        self._ymax = self._vertices[:, 1].max()
        self._zmax = self._vertices[:, 2].max()

        self._xmin = self._vertices[:,0].min()
        self._ymin = self._vertices[:, 1].min()
        self._zmin = self._vertices[:, 2].min()


        self._origin = np.array([x_origin, y_origin, z_origin])

        # TODO: perhaps rewrite this to find the plane

        self._sympy_plane  = Plane(Point3D(self._origin), normal_vector=self._normal)


    def is_intersecting(self,ray):
        """
        checks if ray intersects plane
        
        :param ray: 
        :return: bool, array
        """



        intersecting_point = self._sympy_plane.intersection(ray.sympy_line)[0]

        if 'x' in self._name:

            if self._within_y_bounds(intersecting_point.y) and self._within_z_bounds(intersecting_point.z):
                return True, np.array(map(float, [intersecting_point.x, intersecting_point.y, intersecting_point.z]))



        elif 'y' in self._name:

            if self._within_x_bounds(intersecting_point.x) and self._within_z_bounds(intersecting_point.z):
                return True, np.array(map(float, [intersecting_point.x, intersecting_point.y, intersecting_point.z]))



        elif 'z' in self._name:

            if self._within_y_bounds(intersecting_point.y) and self._within_x_bounds(intersecting_point.x):
                return True, np.array(map(float, [intersecting_point.x, intersecting_point.y, intersecting_point.z]))


        return False, None

    def _within_x_bounds(self,x):

        if x <= self._xmax and self._xmin <= x:

            return True

        else:

            return False


    def _within_y_bounds(self,y):

        if y <= self._ymax and self._ymin <= y:

            return True

        else:

            return False

    def _within_z_bounds(self,z):

        if z <= self._zmax and self._zmin <= z:

            return True

        else:

            return False

    @property
    def origin(self):

        return self._origin



















class Volume(object):
    def __init__(self,name ,x_origin, y_origin, z_origin, height, x_width, y_width, color='grey'):
        self._center = (x_origin, y_origin, z_origin)

        self._color = color

        self._name = name


        self._build_cube(origin=(x_origin - x_width / 2., z_origin - height / 2., y_origin - y_width / 2.),
                         width=x_width,
                         depth=y_width,
                         height=height)

        self._intersections = None




    def _build_cube(self,origin=None, width=1,height=1, depth=1):



        self._planes = collections.OrderedDict()
        for plane in ['+x', '-x', '+y', '-y', '+z', '-z']:

            self._planes[plane] = None



        u, v, w = (0, 0, 0) if origin is None else origin


        grids = []

        for plane in self._planes.keys():


            if '-z' in plane:

                this_grid = self._grid('xy', (u, w), width, depth, v)

                self._planes[plane] = Surface(plane,this_grid[0])

                grids.extend(this_grid)

            if '+z' in plane:

                this_grid = self._grid('xy', (u, w), width, depth, v + height)

                self._planes[plane] = Surface(plane,this_grid[0])

                grids.extend(this_grid)


            if '-y' in plane:


                this_grid = self._grid('xz', (u, v), width, height, w)

                self._planes[plane] = Surface(plane,this_grid[0])

                grids.extend(this_grid)

            if '+y' in plane:

                this_grid = self._grid('xz', (u, v), width, height, w + depth)

                self._planes[plane] = Surface(plane,this_grid[0])

                grids.extend(this_grid)

            if '-x' in plane:

                this_grid = self._grid('yz', (w, v), depth, height, u)

                self._planes[plane] = Surface(plane,this_grid[0])

                grids.extend(this_grid)

            if '+x' in plane:

                this_grid = self._grid('yz', (w, v), depth, height, u + width)

                self._planes[plane] = Surface(plane,this_grid[0])

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
    def planes(self):
        return self._planes

    @property
    def center(self):
        return self._center

    @property
    def name(self):

        return self._name

    def plot(self, ax):
        collection = Poly3DCollection(self._quads,
                                      facecolors=self._color,
                                      alpha=.1)

        collection.set_facecolor((0, 0, 1, .1))

        ax.add_collection3d(collection)

    def intersect_ray(self, ray):

        intersections = collections.OrderedDict()

        for k, v in self._planes.iteritems():

            intersection_info = collections.OrderedDict()

            is_intersecting, point = v.is_intersecting(ray)

            if is_intersecting:

                intersection_info['intersection point'] = point

                # now get the distance between the points

                d2 = (np.power(point - ray.detector_origin, 2)).sum()

                intersection_info['distance'] = np.sqrt(d2)


                intersections[k] = intersection_info


        self._intersections = intersections


    @property
    def all_intersections(self):

        return self._intersections

    @property
    def intersection(self):


        # return the first intersection

        max_distance = 0.

        intersection = None

        for k, v in self._intersections.iteritems():

            if v['distance'] > max_distance:

                intersection = k
                max_distance = v['distance']


        if intersection is None:

            return None, None, None



        return intersection, self._intersections[intersection]['intersection point'], self._intersections[intersection]['distance']






