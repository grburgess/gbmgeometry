import numpy as np


def get_sc_matrix(q1, q2, q3, q4):
    """
    build a sc matrix from quats

    :param q1: 
    :param q2: 
    :param q3: 
    :param q4: 
    :returns: 
    :rtype: 

    """

    sc_matrix = np.zeros((3, 3))

    sc_matrix[0, 0] = q1 ** 2 - q2 ** 2 - q3 ** 2 + q4 ** 2
    sc_matrix[0, 1] = 2.0 * (q1 * q2 + q4 * q3)
    sc_matrix[0, 2] = 2.0 * (q1 * q3 - q4 * q2)
    sc_matrix[1, 0] = 2.0 * (q1 * q2 - q4 * q3)
    sc_matrix[1, 1] = -(q1 ** 2) + q2 ** 2 - q3 ** 2 + q4 ** 2
    sc_matrix[1, 2] = 2.0 * (q2 * q3 + q4 * q1)
    sc_matrix[2, 0] = 2.0 * (q1 * q3 + q4 * q2)
    sc_matrix[2, 1] = 2.0 * (q2 * q3 - q4 * q1)
    sc_matrix[2, 2] = -(q1 ** 2) - q2 ** 2 + q3 ** 2 + q4 ** 2

    return sc_matrix
