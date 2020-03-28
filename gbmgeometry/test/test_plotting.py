from gbmgeometry.utils.plotting.space_plot import animate_in_space, plot_in_space
from gbmgeometry import PositionInterpolator
from gbmgeometry.utils.plotting.sky_point import balrog_to_skypoints
from gbmgeometry.utils.package_utils import get_path_of_data_file


def test_space_plot(interpolator):

    tmin, tmax = interpolator.minmax_time()

    plot_in_space(
        interpolator,
        tmin,
        show_detector_pointing=True,
        show_moon=True,
        show_sun=True,
        show_stars=True,
    )


def test_space_ani(interpolator):

    plot_in_space(
        interpolator,
        10,
        show_detector_pointing=True,
        show_moon=True,
        show_sun=True,
        show_stars=True,
    )


def test_point_space_plotting():

    pi = PositionInterpolator.from_trigdat(get_path_of_data_file("balrog_trig.fits"))

    skypoints = balrog_to_skypoints(
        get_path_of_data_file("balrog.fits"),
        new_nside=2 ** 5,
        cmap="viridis",
        as_point=True,
    )


def test_ray_space_plotting():

    pi = PositionInterpolator.from_trigdat(get_path_of_data_file("balrog_trig.fits"))

    skypoints = balrog_to_skypoints(
        get_path_of_data_file("balrog.fits"),
        new_nside=2 ** 5,
        cmap="viridis",
        as_point=False,
    )

    plot_in_space(
        pi,
        0,
        sky_points=skypoints,
        show_detector_pointing=True,
        show_moon=True,
        earth_time="day",
        show_stars=True,
    )
