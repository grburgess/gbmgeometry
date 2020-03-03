from gbmgeometry.gbm import GBM, get_legal_pairs
from gbmgeometry.gbm_detector import BGO0, BGO1
from gbmgeometry.gbm_detector import NaI0, NaI1, NaI2, NaI3, NaI4, NaI5
from gbmgeometry.gbm_detector import NaI6, NaI7, NaI8, NaI9, NaIA, NaIB
from gbmgeometry.gbm_frame import GBMFrame
from gbmgeometry.position_interpolator import PositionInterpolator
from gbmgeometry.spacecraft.fermi import Fermi


gbm_detector_list = {
    "n0": NaI0,
    "n1": NaI1,
    "n2": NaI2,
    "n3": NaI3,
    "n4": NaI4,
    "n5": NaI5,
    "n6": NaI6,
    "n7": NaI7,
    "n8": NaI8,
    "n9": NaI9,
    "na": NaIA,
    "nb": NaIB,
    "b0": BGO0,
    "b1": BGO1,
}


from gbmgeometry.utils.gbm_time import GBMTime

__all__ = [
    "GBMFrame",
    "GBM",
    "NaI0",
    "NaI1",
    "NaI2",
    "NaI3",
    "NaI4",
    "NaI5",
    "NaI6",
    "NaI7",
    "NaI8",
    "NaI9",
    "NaIA",
    "NaIB",
    "BGO0",
    "BGO1",
    "PositionInterpolator",
    "get_legal_pairs",
    "Fermi",
    "GBMTime",
    "gbm_detector_list",
]
