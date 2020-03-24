import numpy as np
import ipyvolume as ipv
from gbmgeometry.utils.plotting.heavenly_bodies import (
    Sol,
    Moon,
    Earth,
    Sphere,
    StarField,
)
from gbmgeometry.gbm import GBM


def compute_distance(x, y, z, radius):

    dist = np.sqrt(x * x + y * y + z * z)

    dist += radius

    return dist


class FermiPoint(Sphere):
    def __init__(self, x, y, z, color="#10DC9B"):
        """
        A dummy point for fermi to keep sizes in check
        """

        super(FermiPoint, self).__init__(
            ax=None, x=x, y=y, z=z, detail_level=20, radius=100.0, color=color,
        )


def animate_in_space(
    position_interpolator,
    n_step=200,
    show_detector_pointing=False,
    show_earth=True,
    show_sun=False,
    show_moon=False,
    background_color="#070323",
    detector_scaling_factor=20000.0,
    show_stars=False,
):
    """
    Animiate fermi in Space!

    :param position_interpolator: 
    :param n_step: 
    :param show_detector_pointing: 
    :param show_earth: 
    :param show_sun: 
    :param show_moon: 
    :param background_color: 
    :param detector_scaling_factor: 
    :param show_stars: 
    :returns: 
    :rtype: 

    """

    fig = ipv.figure()

    ipv.pylab.style.box_off()
    ipv.pylab.style.axes_off()
    ipv.pylab.style.set_style_dark()
    ipv.pylab.style.background_color(background_color)

    tmin, tmax = position_interpolator.minmax_time()

    time = np.linspace(tmin, tmax, n_step)

    artists = []

    distances = [8000]

    if show_earth:

        earth = Earth()

        earth.plot()

    if show_sun:

        xs = []
        ys = []
        zs = []

        for t in time:

            sun_pos = position_interpolator.sun_position(t)
            x, y, z = sun_pos.cartesian.xyz.to("km").value

            xs.append(x)
            ys.append(y)
            zs.append(z)

        sol = Sol(np.array(xs), np.array(ys), np.array(zs))

        distances.append(compute_distance(x, y, z, sol.radius))

        artists.append(sol.plot())

    if show_moon:

        xs = []
        ys = []
        zs = []

        for t in time:

            moon_pos = position_interpolator.moon_position(t)
            x, y, z = moon_pos.cartesian.xyz.to("km").value

            xs.append(x)
            ys.append(y)
            zs.append(z)

        moon = Moon(np.array(xs), np.array(ys), np.array(zs))
        distances.append(compute_distance(x, y, z, moon.radius))
        artists.append(moon.plot())

    # now get fermi position
    sxs = []
    sys = []
    szs = []

    if show_detector_pointing:

        distances.append(detector_scaling_factor)

        gbm = GBM(
            position_interpolator.quaternion(tmin), position_interpolator.sc_pos(tmin),
        )

        dets_x = {}
        dets_y = {}
        dets_z = {}

        for k, _ in gbm.detectors.items():

            dets_x[k] = []
            dets_y[k] = []
            dets_z[k] = []

    for t in time:

        sx, sy, sz = position_interpolator.sc_pos(t)

        sxs.append(sx)
        sys.append(sy)
        szs.append(sz)

        if show_detector_pointing:

            gbm.set_quaternion(position_interpolator.quaternion(t))

            for k, v in gbm.detectors.items():

                x, y, z = v.center_icrs.cartesian.xyz.value * max(distances)

                dets_x[k].append([sx, x])
                dets_y[k].append([sy, y])
                dets_z[k].append([sz, z])

    if show_detector_pointing:

        for k, v in gbm.detectors.items():

            dets_x[k] = np.array(dets_x[k])
            dets_y[k] = np.array(dets_y[k])
            dets_z[k] = np.array(dets_z[k])

            if k.startswith("b"):
                color = "#EE3A00"

            else:
                color = "#EEEE01"

            artists.append(ipv.pylab.plot(dets_x[k], dets_y[k], dets_z[k], color=color))

    sxs = np.array(sxs)
    sys = np.array(sys)
    szs = np.array(szs)

    fermi = FermiPoint(sxs, sys, szs)
    artists.append(fermi.plot())

    if show_stars:

        sf = StarField(n_stars=200, distance=max(distances) - 2)
        sf.plot()

    ipv.xyzlim(max(distances))

    ipv.animation_control(artists)

    ipv.show()


def plot_in_space(
    position_interpolator,
    time,
    show_detector_pointing=False,
    show_earth=True,
    show_sun=False,
    show_moon=False,
    background_color="#070323",
    detector_scaling_factor=20000.0,
    show_stars=False,
):
    """
    Plot Fermi in Space!

    :param position_interpolator: 
    :param time: 
    :param show_detector_pointing: 
    :param show_earth: 
    :param show_sun: 
    :param show_moon: 
    :param background_color: 
    :param detector_scaling_factor: 
    :returns: 
    :rtype: 

    """

    fig = ipv.figure()

    ipv.pylab.style.box_off()
    ipv.pylab.style.axes_off()
    ipv.pylab.style.set_style_dark()
    ipv.pylab.style.background_color(background_color)

    distances = [8000]

    if show_earth:

        earth = Earth()

        earth.plot()

    if show_sun:

        sun_pos = position_interpolator.sun_position(time)
        x, y, z = sun_pos.cartesian.xyz.to("km").value

        sol = Sol(x, y, z)
        distances.append(compute_distance(x, y, z, sol.radius))
        sol.plot()

    if show_moon:

        moon_pos = position_interpolator.moon_position(time)
        x, y, z = moon_pos.cartesian.xyz.to("km").value

        moon = Moon(x, y, z)
        distances.append(compute_distance(x, y, z, moon.radius))
        moon.plot()

    # now get fermi position

    sx, sy, sz = position_interpolator.sc_pos(time)

    fermi = FermiPoint(sx, sy, sz)
    fermi.plot()

    if show_detector_pointing:

        distances.append(detector_scaling_factor)

        gbm = GBM(
            position_interpolator.quaternion(time), position_interpolator.sc_pos(time),
        )

        for k, v in gbm.detectors.items():

            x, y, z = v.center_icrs.cartesian.xyz.value * max(distances)

            x_line = np.array([sx, x])
            y_line = np.array([sy, y])
            z_line = np.array([sz, z])

            if k.startswith("b"):
                color = "#EE3A00"

            else:
                color = "#EEEE01"

            ipv.pylab.plot(x_line, y_line, z_line, color=color)

    if show_stars:

        sf = StarField(n_stars=100, distance=max(distances) - 2)
        sf.plot()

    ipv.xyzlim(max(distances))

    ipv.show()
