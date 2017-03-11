import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import sys

# Bad fits throw exceptions
np.seterr(all='raise')


def hyperbola(concentrations, f0, delta_f, kd):
    """
    :param concentrations: array of titrant concentrations
    :param f0: Y-intercept
    :param delta_f: Maximum change in signal
    :param kd: Dissociation constant
    :return: array of Y values

    """
    return f0 + ((delta_f * concentrations) / (concentrations + kd))


def fit_hyperbola(concentrations, signals):
    if len(concentrations) < 3 or len(signals) < 3:
        raise ValueError("You must have at least two measurements for a fit.")
    if len(set(concentrations)) != len(concentrations):
        raise ValueError("You may not provide multiple measurements for the same X-axis value.")
    if len(concentrations) != len(signals):
        raise ValueError("You must have a signal for each given concentration!")
    if np.any(np.isnan(signals)) or np.any(np.isnan(concentrations)):
        raise ValueError("Your data may not contain NaNs")
    if np.any(concentrations[concentrations < 10**-50]) or np.any(concentrations[concentrations > 10**50]):
        raise ValueError("You provided unreasonably small or large values.")
    if np.max(signals) - np.min(signals) < 10**-50:
        raise ValueError("Your signal values do not change sufficiently for a fit.")
    try:
        slope, intercept, _, _, _ = stats.linregress(concentrations, signals)
    except FloatingPointError:
        raise ValueError("Your signals make no sense.")
    f0_min, f0_max = (0.0, np.inf) if intercept > 0 else (-np.inf, 0.0)
    delta_f_min, delta_f_max = (0.0, np.inf) if slope > 0 else (-np.inf, 0.0)
    (f0, delta_f, kd), covariance = curve_fit(hyperbola, concentrations, signals, bounds=((f0_min, delta_f_min, 0.0),
                                                                                          (f0_max, delta_f_max, np.inf)))
    f0_stddev = covariance[0, 0] ** 0.5
    delta_f_stddev = covariance[1, 1] ** 0.5
    kd_stddev = covariance[2, 2] ** 0.5
    return f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev


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


def fit_quadratic(concentrations, signals):
    slope, intercept, _, _, _ = stats.linregress(concentrations, signals)
    f0_min, f0_max = (0.0, np.inf) if intercept > 0 else (-np.inf, 0.0)
    delta_f_min, delta_f_max = (0.0, np.inf) if slope > 0 else (-np.inf, 0.0)
    (f0, delta_f, kd, constant), covariance = curve_fit(quadratic, concentrations, signals, bounds=((f0_min, delta_f_min, 0.0, 0.0),
                                                                                                    (f0_max, delta_f_max, np.inf, np.inf)))
    f0_stddev = covariance[0, 0] ** 0.5
    delta_f_stddev = covariance[1, 1] ** 0.5
    kd_stddev = covariance[2, 2] ** 0.5
    constant_stddev = covariance[3,3] ** 0.5
    return f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev, constant, constant_stddev
