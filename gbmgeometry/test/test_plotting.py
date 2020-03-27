from gbmgeometry.utils.plotting.space_plot import animate_in_space, plot_in_space


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
