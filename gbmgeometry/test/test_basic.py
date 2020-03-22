from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u
from gbmgeometry import PositionInterpolator, NaIA, GBMFrame, GBM
from gbmgeometry.utils.package_utils import get_path_of_data_file
from gbmgeometry.utils.plotting.skyplot import skyplot
import numpy as np



def test_interp():

    trigdat = get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit")

    interp_trig = PositionInterpolator.from_trigdat(trigdat)

    interp_trig.quaternion(0)

    interp_trig.sc_pos(0)

    trigdat_h5 = get_path_of_data_file("trigdat.h5")

    interp_trig_h5 = PositionInterpolator.from_trigdat_hdf5(trigdat_h5)

    assert np.all(interp_trig_h5.quaternion(0) == interp_trig.quaternion(0))

    assert np.all(interp_trig_h5.sc_pos(0) == interp_trig.sc_pos(0))

    poshist = get_path_of_data_file("glg_poshist_all_151013_v00.fit")

    interp_pos = PositionInterpolator.from_poshist(poshist)

    interp_pos.quaternion(interp_pos.time[0])

    interp_pos.sc_pos(interp_pos.time[0])

    poshist_h5 = get_path_of_data_file("posthist.h5")

    interp_pos_h5 = PositionInterpolator.from_poshist_hdf5(poshist_h5)

    assert np.all(interp_pos_h5.quaternion(interp_pos.time_h5[0]) == interp_pos.quaternion(interp_pos.time[0]))

    assert np.all(interp_pos_h5.sc_pos(interp_pos_h5.time[0]) == interp_pos.sc_pos(interp_pos.time[0]))


def test_detector(na, interpolator):

    na.plot_pointing(fov=5)

    na.get_center()

    na.set_quaternion(interpolator.quaternion(100))

    na.get_center()

    na.set_quaternion(interpolator.quaternion(100))


def test_coord_change(na, interpolator):


    center_j2000 = na.get_center().icrs

    q1, q2, q3, q4 = interpolator.quaternion(100.0)

    center_j2000.transform_to(
        GBMFrame(quaternion_1=q1, quaternion_2=q2, quaternion_3=q3, quaternion_4=q4)
    )


    center_j2000 = na.get_center().icrs

    out = interpolator.quaternion_dict(100.0)

    center_j2000.transform_to(
        GBMFrame(**out)
    )


    

def test_all_gbm(interpolator):

    myGBM = GBM(interpolator.quaternion(0), sc_pos=interpolator.sc_pos(0) * u.km)

    myGBM.get_centers()

    [x.icrs for x in myGBM.get_centers()]

    ax = skyplot(background_color="#47496C")

    for k, v in myGBM.detectors.items():

        v.plot_pointing(
            ax=ax, fov=5, edgecolor="limegreen", lw=1, facecolor="limegreen", alpha=0.4
        )
