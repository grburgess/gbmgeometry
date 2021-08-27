import os
import re
from urllib.error import HTTPError

import astropy.io.fits as fits

from gbmgeometry.utils.file_converters import (convert_poshist2hdf5,
                                               convert_trigdat2hdf5)

# base url

_base_url = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm"
_trigger_base_url = "/".join([_base_url, "triggers"])
_daily_base_url = "/".join([_base_url, "daily"])

burst_number_match = re.compile("^(bn|GRB)?(\d{2})(\d{2})(\d{2})(\d{3})")


def download_trigdat(burst_number, version=None, destination="."):
    """
    download the trigdat data for a given GRB trigger as an
    HDF5 file

    :param burst_number: grb|bn YYMMDDXXX
    :param version:
    :param destination: where this file should go
    :returns:
    :rtype:

    """
    # use regex to find the needed info from the burst number

    groups = burst_number_match.match(burst_number).groups()

    if len(groups) == 5:
        _, year, month, day, frac = groups

    elif len(groups) == 4:

        year, month, day, frac = groups

    else:

        raise RuntimeError(
            f"{burst_number} is not in the correct format! Should look like bn200407388"
        )

    full_year = f"20{year}"
    proper_bn = f"bn{year}{month}{day}{frac}"

    # if a version is specified

    if version is not None:

        versions = [version]

    # else we will try until success

    else:

        versions = range(4)

    found = False

    for version in versions:

        # assemble the URL

        url = "/".join(
            [
                _trigger_base_url,
                full_year,
                proper_bn,
                "current",
                f"glg_trigdat_all_{proper_bn}_v0{int(version)}.fit",
            ]
        )

        out_file = os.path.join(
            destination, f"trigdat_{proper_bn}_v0{int(version)}.h5")

        try:

            convert_trigdat2hdf5(url, out_file)

            found = True

            break

        except (HTTPError):

            pass

    if not found:

        print("Sorry there were no trigdat data for this GRB")

    else:

        print(f"The data have been downloaded to {out_file}")

    return out_file


def download_posthist(year, month, day, destination="."):
    """
    download the position history file from GBM for the
    given day as an HDF5 file

    :param year: '20'
    :param month: '05'
    :param day: '12'
    :param destination: where this file should go
    :returns:
    :rtype:

    """

    full_year = f"20{year}"

    # assemble the URL

    url = "/".join(
        [
            _daily_base_url,
            full_year,
            month,
            day,
            "current",
            f"glg_poshist_all_{year}{month}{day}_v00.fit",
        ]
    )

    
    out_file = os.path.join(destination, f"poshist_{year}{month}{day}.h5")

    try:

        convert_poshist2hdf5(url, out_file)

        found = True

    except (HTTPError):


        try:

            url = "/".join(
                [
                    _daily_base_url,
                    full_year,
                    month,
                    day,
                    "current",
                    f"glg_poshist_all_{year}{month}{day}_v01.fit",
                ]
            )
            
            convert_poshist2hdf5(url, out_file)

            found = True

        except(HTTPError):

            print("Sorry there were no posthist data for this the day")

            return
            
    print(f"The data have been downloaded to {out_file}")

    return out_file


def get_official_location(burst_number, version=None):
    """
    download the trigdat data for a given GRB trigger as an
    HDF5 file

    :param burst_number: grb|bn YYMMDDXXX
    :param version:
    :param destination: where this file should go
    :returns:
    :rtype:

    """
    # use regex to find the needed info from the burst number

    groups = burst_number_match.match(burst_number).groups()

    if len(groups) == 5:
        _, year, month, day, frac = groups

    elif len(groups) == 4:

        year, month, day, frac = groups

    else:

        raise RuntimeError(
            f"{burst_number} is not in the correct format! Should look like bn200407388"
        )

    full_year = f"20{year}"
    proper_bn = f"bn{year}{month}{day}{frac}"

    # if a version is specified

    if version is not None:

        versions = [version]

    # else we will try until success

    else:

        versions = range(20)

 

    for version in versions[::-1]:

        # assemble the URL

        url = "/".join(
            [
                _trigger_base_url,
                full_year,
                proper_bn,
                "current",
                f"glg_tcat_all_{proper_bn}_v{str(version).zfill(2)}.fit",
            ]
        )

        try:

            with fits.open(url) as f:

                header = f[0].header

            break

        except:

            pass

    return (header["RA_OBJ"], header["DEC_OBJ"], header["ERR_RAD"])
