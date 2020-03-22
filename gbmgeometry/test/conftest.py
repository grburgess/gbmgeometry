from gbmgeometry import PositionInterpolator, NaIA, GBMFrame, GBM
from gbmgeometry.utils.package_utils import get_path_of_data_file

import pytest


@pytest.fixture(scope="session")
def interpolator():

    trigdat_h5 = get_path_of_data_file("trigdat.h5")

    interp_trig_h5 = PositionInterpolator.from_trigdat_hdf5(trigdat_h5)

    return interp_trig_h5


@pytest.fixture(scope="session")
def na(interpolator):

    na = NaIA(interpolator.quaternion(1), sc_pos=interpolator.sc_pos(1))

    return na
