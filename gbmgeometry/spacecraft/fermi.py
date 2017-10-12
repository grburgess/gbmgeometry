from lat import LAT, LATRadiatorMinus, LATRadiatorPlus
from solar_panels import SolarPanelMinus, SolarPanelPlus
from gbmgeometry.gbm import GBM

from geometry import Ray

import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np



class Fermi(object):

    def __init__(self, quaternion, sc_pos=None):

        # build fermi

        self._lat = LAT()

        self._lat_radiator_plus = LATRadiatorPlus()
        self._lat_radiator_minus = LATRadiatorMinus()

        self._solar_panel_plus = SolarPanelPlus()
        self._solar_panel_minus = SolarPanelMinus()

        self._gbm = GBM(quaternion,sc_pos)


        # attach the components to fermi

        self._spacecraft_components = collections.OrderedDict()

        self._spacecraft_components[self._lat.name] = self._lat

        self._spacecraft_components[self._lat_radiator_minus.name] = self._lat_radiator_minus
        self._spacecraft_components[self._lat_radiator_plus.name] = self._lat_radiator_plus

        self._spacecraft_components[self._solar_panel_plus.name] = self._solar_panel_plus
        self._spacecraft_components[self._solar_panel_minus.name] = self._solar_panel_minus


        # add lists for each detector to rays

        self._rays = collections.OrderedDict()

        for name in self._gbm.detectors.iterkeys():

            self._rays[name] = []


    @property
    def spacecraft_components(self):

        return self._spacecraft_components


    @property
    def rays(self):

        return self._rays


    def add_ray(self, ray_coordinate,probability=None):


        for name, det in self._gbm.detectors.iteritems():


            ray = Ray(det, ray_coordinate, probability)

            self._rays[name].append(ray)




    def compute_intersections(self):




    def plot_fermi(self, ax=None, with_rays=False):


        if ax is None:

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        else:

            fig = ax.get_figure()




        for name, component in self._spacecraft_components.iteritems():

            component.plot(ax)


        for name, det in self._gbm.detectors.iteritems():

            ax.scatter(*det.mount_point,color='#FFC300')
            ax.text3D(*det.mount_point,s=name)



        if with_rays:

            for det in self._rays.itervalues():


                for ray in det:

                    ray.plot(ax)




        ax.set_xlabel('SCX')
        ax.set_ylabel('SCY')
        ax.set_zlabel('SCZ')

        ax.set_zlim(0, 300)
        ax.set_xlim(-400, 400)
        ax.set_ylim(-400, 400)

        ax.grid(False)
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')

        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False





        return fig










