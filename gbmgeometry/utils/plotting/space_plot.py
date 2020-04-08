import numpy as np
import ipyvolume as ipv
from gbmgeometry.utils.plotting.heavenly_bodies import (
    Sol,
    Moon,
    Earth,
    StarField,
)
from gbmgeometry.gbm import GBM
from gbmgeometry.spacecraft.fermi import Fermi
from gbmgeometry.geometry import Sphere
from gbmgeometry.geometry.cone import Cone


_open_angle = 10.0

_det_colors = dict(
    n0="#CC3311",
    n1="#CC3311",
    n2="#CC3311",
    n3="#009988",  # teal
    n4="#009988",
    n5="#009988",
    n6="#EE7733",
    n7="#EE7733",
    n8="#EE7733",
    n9="#0077BB",
    na="#0077BB",
    nb="#0077BB",
    b0="#F2E300",
    b1="#F2E300",
)


def compute_distance(x, y, z, radius):

    dist = np.sqrt(x * x + y * y + z * z)

    dist += radius

    return dist


def animate_in_space(
    position_interpolator,
    n_step=200,
    show_detector_pointing=False,
    show_earth=True,
    show_sun=False,
    show_moon=False,
    background_color="#01000F",
    detector_scaling_factor=20000.0,
    show_stars=False,
    show_inactive=False,
    earth_time="night",
    realistic=True,
    interval=200,
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

    distances = [15000]

    if show_earth:

        astro_times = [position_interpolator.astro_time(t) for t in time]

        earth = Earth(
            earth_time=earth_time, realistic=realistic, astro_time=astro_times,
        )

        tmp = earth.plot()

        if realistic:
            artists.append(tmp)

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

    x_off = []
    y_off = []
    z_off = []

    if show_detector_pointing:

        distances.append(detector_scaling_factor)

        gbm = GBM(
            position_interpolator.quaternion(tmin), position_interpolator.sc_pos(tmin),
        )

        dets_xo = {}
        dets_yo = {}
        dets_zo = {}
        dets_xp = {}
        dets_yp = {}
        dets_zp = {}

        for k, _ in gbm.detectors.items():

            dets_xo[k] = []
            dets_yo[k] = []
            dets_zo[k] = []
            dets_xp[k] = []
            dets_yp[k] = []
            dets_zp[k] = []

    for t in time:

        sx, sy, sz = position_interpolator.sc_pos(t)

        sxs.append(sx)
        sys.append(sy)
        szs.append(sz)

        if not position_interpolator.is_fermi_active(t):
            x_off.append(sx)
            y_off.append(sy)
            z_off.append(sz)

        if show_detector_pointing:

            gbm.set_quaternion(position_interpolator.quaternion(t))

            for k, v in gbm.detectors.items():

                x, y, z = v.center_icrs.cartesian.xyz.value * max(distances)

                dets_xo[k].append(sx)
                dets_yo[k].append(sy)
                dets_zo[k].append(sz)
                dets_xp[k].append(sx + x)
                dets_yp[k].append(sy + y)
                dets_zp[k].append(sz + z)

    if show_detector_pointing:

        for k, v in gbm.detectors.items():

            dets_xo[k] = np.array(dets_xo[k])
            dets_yo[k] = np.array(dets_yo[k])
            dets_zo[k] = np.array(dets_zo[k])
            dets_xp[k] = np.array(dets_xp[k])
            dets_yp[k] = np.array(dets_yp[k])
            dets_zp[k] = np.array(dets_zp[k])

            color = _det_colors[k]

            cone = Cone(
                dets_xo[k],
                dets_yo[k],
                dets_zo[k],
                dets_xp[k],
                dets_yp[k],
                dets_zp[k],
                _open_angle,
            )
            artists.append(cone.plot(color=color))

    #            artists.append(ipv.pylab.plot(dets_x[k], dets_y[k], dets_z[k], color=color))

    sxs = np.array(sxs)
    sys = np.array(sys)
    szs = np.array(szs)
    if show_inactive:

        ipv.pylab.scatter(
            np.array(x_off),
            np.array(y_off),
            np.array(z_off),
            color="#DC1212",
            alpha=0.5,
            marke="circle_2d",
            size=1,
        )

    # fermi = FermiPoint(sxs, sys, szs)
    # artists.append(fermi.plot())

    fermi_real = Fermi(
        position_interpolator.quaternion(time),
        sc_pos=position_interpolator.sc_pos(time),
        transform_to_space=True,
    )
    artists.extend(fermi_real.plot_fermi_ipy())

    if show_stars:

        sf = StarField(n_stars=200, distance=max(distances) - 2)
        sf.plot()

    ipv.xyzlim(max(distances))

    ipv.animation_control(artists, interval=interval)

    ipv.show()

    return fig


def plot_in_space(
    position_interpolator,
    time,
    show_detector_pointing=False,
    show_earth=True,
    show_sun=False,
    show_moon=False,
    background_color="#01000F",
    detector_scaling_factor=20000.0,
    show_stars=False,
    show_orbit=True,
    realistic=True,
    earth_time="night",
    sky_points=None,
    min_distance=8000,
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

    fig = ipv.figure(width=800, height=600)

    ipv.pylab.style.box_off()
    ipv.pylab.style.axes_off()
    ipv.pylab.style.set_style_dark()
    ipv.pylab.style.background_color(background_color)

    distances = [min_distance]

    if sky_points is not None:
        sky_points = np.atleast_1d(sky_points)

    if show_orbit:

        tmin, tmax = position_interpolator.minmax_time()
        tt = np.linspace(tmin, tmax, 500)

        sc_pos = position_interpolator.sc_pos(tt)

        ipv.plot(sc_pos[:, 0], sc_pos[:, 1], sc_pos[:, 2], lw=0.5)

    if show_earth:

        earth = Earth(
            earth_time=earth_time,
            realistic=realistic,
            astro_time=position_interpolator.astro_time(time),
        )

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

        moon = Moon(x, y, z, realistic=True)
        distances.append(compute_distance(x, y, z, moon.radius))
        moon.plot()

    # now get fermi position

    sx, sy, sz = position_interpolator.sc_pos(time)

    fermi_real = Fermi(
        position_interpolator.quaternion(time),
        sc_pos=position_interpolator.sc_pos(time),
        transform_to_space=True,
    )
    fermi_real.plot_fermi_ipy()

    if show_detector_pointing:

        distances.append(detector_scaling_factor)

        gbm = GBM(
            position_interpolator.quaternion(time), position_interpolator.sc_pos(time),
        )

        for k, v in gbm.detectors.items():
            x, y, z = v.center_icrs.cartesian.xyz.value * max(distances)

            # x_line = np.array([sx, sx + x])
            # y_line = np.array([sy, sy + y])
            # z_line = np.array([sz, sz + z])

            color = _det_colors[k]

            cone = Cone(sx, sy, sz, sx + x, sy + y, sz + z, _open_angle)
            cone.plot(color=color)

            # ipv.pylab.plot(x_line, y_line, z_line, color=color)

    if sky_points is not None:
        for sp in sky_points:

            sp.plot(sx, sy, sz, max(distances))

            # distances.append(sp.distance)

    if show_stars:

        sf = StarField(n_stars=100, distance=max(distances) - 2)
        sf.plot()

    # move the camera here

    fig.camera.position = tuple(
        position_interpolator.sc_pos(time)
        / np.linalg.norm(position_interpolator.sc_pos(time))
    )

    ipv.xyzlim(max(distances))

    ipv.show()

    return fig
