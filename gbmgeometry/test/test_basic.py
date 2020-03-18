from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u
from gbmgeometry import *
from gbmgeometry.utils.package_utils import get_path_of_data_file
from gbmgeometry.utils.plotting.skyplot import skyplot

def test_interp():

    trigdat = get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit")

    interp = PositionInterpolator(trigdat=trigdat)

    interp.quaternion(0)

    interp.sc_pos(0)


def test_detector():

    trigdat = get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit")

    interp = PositionInterpolator(trigdat=trigdat)

    na = NaIA(interp.quaternion(1))

    na.plot_pointing(fov=5)
    
    na.get_center()

    na.set_quaternion(interp.quaternion(100))

    na = NaIA(interp.quaternion(1), sc_pos=interp.sc_pos(1))

    na.get_center()

    na.set_quaternion(interp.quaternion(100))


def test_coord_change():

    trigdat = get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit")

    interp = PositionInterpolator(trigdat=trigdat)

    na = NaIA(interp.quaternion(1))

    center_j2000 = na.get_center().icrs

    q1, q2, q3, q4 = interp.quaternion(100.0)

    center_j2000.transform_to(
        GBMFrame(quaternion_1=q1, quaternion_2=q2, quaternion_3=q3, quaternion_4=q4)
    )


def test_all_gbm():

    trigdat = get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit")

    interp = PositionInterpolator(trigdat=trigdat)

    myGBM = GBM(interp.quaternion(0), sc_pos=interp.sc_pos(0) * u.km)

    myGBM.get_centers()

    [x.icrs for x in myGBM.get_centers()]

    ax = skyplot(background_color="#47496C")

    for k,v in myGBM.detectors.items():
    
        v.plot_pointing(ax=ax,fov=5, edgecolor='limegreen', lw=1, facecolor='limegreen', alpha=.4)

