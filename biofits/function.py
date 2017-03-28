"""
Mathematical functions we might wish to fit our data to.

"""
import numpy as np


def hyperbola(concentrations, yint, delta_y, kd):
    """
    :param concentrations: array of titrant concentrations
    :param yint: Y-intercept
    :param delta_y: Total change in signal
    :param kd: Dissociation constant
    :return: array of Y values

    """
    return yint + ((delta_y * concentrations) / (concentrations + kd))


def quadratic(concentrations, yint, delta_y, kd, constant):
    """
    :param concentrations: array of titrant concentrations
    :param yint: Y-intercept
    :param delta_y: Total change in signal
    :param kd: Dissociation constant
    :param constant: Concentration of substance that is held constant during a titration

    :return: array of Y values

    """
    b = constant + concentrations + kd
    return yint + delta_y * ((b - np.sqrt(np.power(b, 2) - 4 * constant * concentrations)) / (2 * constant))
