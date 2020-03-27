import h5py
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


def convert_poshist2hdf5(src_file, dest_file):
    """
    convert a GBM poshist file to HDF5

    :param src_file: poshist file
    :param dest_file: destination file
    :returns: 
    :rtype: 

    """

    with fits.open(src_file) as poshist:

        time = poshist["GLAST POS HIST"].data["SCLK_UTC"]

        flags = poshist["GLAST POS HIST"].data["FLAGS"]

        quats = np.array(
            [
                poshist["GLAST POS HIST"].data["QSJ_1"],
                poshist["GLAST POS HIST"].data["QSJ_2"],
                poshist["GLAST POS HIST"].data["QSJ_3"],
                poshist["GLAST POS HIST"].data["QSJ_4"],
            ]
        ).T

        sc_pos = np.array(
            [
                poshist["GLAST POS HIST"].data["POS_X"],
                poshist["GLAST POS HIST"].data["POS_Y"],
                poshist["GLAST POS HIST"].data["POS_Z"],
            ]
        ).T

    with h5py.File(dest_file, "w") as f:

        f.create_dataset("time", data=time, compression="lzf")
        f.create_dataset("quats", data=quats, compression="lzf")
        f.create_dataset("sc_pos", data=sc_pos, compression="lzf")
        f.create_dataset("flags", data=flags, compression="lzf")


def convert_trigdat2hdf5(src_file, dest_file):
    """
    convert a GBM trigdat file to HDF5

    :param src_file: poshist file
    :param dest_file: destination file
    :returns: 
    :rtype: 
    """

    with fits.open(src_file) as trigdat:

        trigtime = trigdat["EVNTRATE"].header["TRIGTIME"]
        tstart = trigdat["EVNTRATE"].data["TIME"] - trigtime

        trigtime = trigtime

        quats = trigdat["EVNTRATE"].data["SCATTITD"]
        sc_pos = trigdat["EVNTRATE"].data["EIC"]

        sort_mask = np.argsort(tstart)
        tstart = tstart[sort_mask]

        quats = quats[sort_mask]
        sc_pos = sc_pos[sort_mask]

        time = tstart

    with h5py.File(dest_file, "w") as f:

        f.create_dataset("time", data=time, compression="lzf")
        f.create_dataset("quats", data=quats, compression="lzf")
        f.create_dataset("sc_pos", data=sc_pos, compression="lzf")

        f.attrs["trigtime"] = trigtime


def make_image_hdf5(image, file_name, detail_level=5):
    img = plt.imread(image)
    # define a grid matching the map size, subsample along with pixels
    theta = np.linspace(0, np.pi, img.shape[0])
    phi = np.linspace(0, 2 * np.pi, img.shape[1])

    count = int(detail_level * 180)
    theta_inds = np.linspace(0, img.shape[0] - 1, count).round().astype(int)
    phi_inds = np.linspace(0, img.shape[1] - 1, count).round().astype(int)
    theta = theta[theta_inds]
    phi = phi[phi_inds]
    img = img[np.ix_(theta_inds, phi_inds)]

    theta, phi = np.meshgrid(theta, phi)

    # sphere
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    color = img / 255

    with h5py.File(file_name, "w") as f:

        f.create_dataset("x", data=x, compression="lzf")
        f.create_dataset("y", data=y, compression="lzf")
        f.create_dataset("z", data=z, compression="lzf")
        f.create_dataset("color", data=color, compression="lzf")
