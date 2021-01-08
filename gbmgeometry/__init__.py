from gbmgeometry.gbm import GBM  # , get_legal_pairs
from gbmgeometry.gbm_detector import (BGO0, BGO1, NaI0, NaI1, NaI2, NaI3, NaI4,
                                      NaI5, NaI6, NaI7, NaI8, NaI9, NaIA, NaIB)
from gbmgeometry.gbm_frame import GBMFrame
from gbmgeometry.io.data_download import (download_posthist, download_trigdat,
                                          get_official_location)
from gbmgeometry.position_interpolator import PositionInterpolator
from gbmgeometry.spacecraft.fermi import Fermi
from gbmgeometry.utils.file_converters import (convert_poshist2hdf5,
                                               convert_trigdat2hdf5)
from gbmgeometry.utils.gbm_time import GBMTime
from gbmgeometry.utils.plotting.space_plot import (animate_in_space,
                                                   plot_in_space)

from ._version import get_versions

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


__all__ = [
    "plot_in_space",
    "animate_in_space",
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
    "convert_poshist2hdf5",
    "convert_trigdat2hdf5",
    "PositionInterpolator",
    #    "get_legal_pairs",
    "Fermi",
    "GBMTime",
    "gbm_detector_list",
]


__version__ = get_versions()["version"]
del get_versions
