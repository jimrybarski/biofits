import numpy as np
from scipy.optimize import curve_fit


def hyperbola(concentrations, f0, delta_f, kd):
    """
    :param concentrations: array of titrant concentrations
    :param f0: Y-intercept
    :param delta_f: Maximum change in signal
    :param kd: Dissociation constant
    :return: array of Y values

    """
    return f0 + ((delta_f * concentrations) / (concentrations + kd))


def quadratic(concentrations, f0, delta_f, kd, constant):
    """
    :param concentrations: array of titrant concentrations
    :param f0: Y-intercept
    :param delta_f: Maximum change in signal
    :param kd: Dissociation constant
    :param constant: Concentration of molecule that is held constant during a titration
    :return: array of Y values

    """
    b = constant + concentrations + kd
    return f0 + delta_f * ((b - np.sqrt(np.power(b, 2) - 4 * constant * concentrations)) / (2 * constant))


def fit_hyperbola(concentrations, signals):
    (f0, delta_f, kd), pcov = curve_fit(hyperbola, concentrations, signals)
    f0_stddev = pcov[0, 0]**0.5
    delta_f_stddev = pcov[1, 1]**0.5
    kd_stddev = pcov[2, 2]**0.5


def fit_quadratic(concentrations, signals):
    (f0, delta_f, kd), pcov = curve_fit(quadratic, concentrations, signals)
    f0_stddev = pcov[0, 0] ** 0.5
    delta_f_stddev = pcov[1, 1] ** 0.5
    kd_stddev = pcov[2, 2] ** 0.5
