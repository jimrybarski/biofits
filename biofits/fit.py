import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
from biofits.function import hyperbola, quadratic, exponential
import sys

# Bad fits throw exceptions
np.seterr(all='raise')


def fit_hyperbola(concentrations, signals):
    """
    :param concentrations: X-axis values representing concentrations in arbitrary units
    :param signals: Y-axis values representing some kind of signal. Don't normalize this.

    :return:
        yint: the Y-intercept of the fit, often the background signal
        yint_stddev: standard deviation of the error of yint
        delta_y: total height of the fit
        delta_y_stddev: standard deviation of the error of delta_y
        kd: the dissociation constant
        kd_stddev: the standard deviation of the error of kd

    """
    slope, intercept, _, _, _ = stats.linregress(concentrations, signals)
    yint_min, yint_max = (0.0, np.inf) if intercept > 0 else (-np.inf, 0.0)
    delta_y_min, delta_y_max = (0.0, np.inf) if slope > 0 else (-np.inf, 0.0)

    sigmas = [np.std(s) for s in signals]
    if 0 in sigmas:
        # if we have scalar values, the standard deviation will be zero everywhere and this will cause a division by
        # zero to occur. In this case, we just weight all points equally
        sigmas = np.ones((len(signals),))
    (yint, delta_y, kd), covariance = curve_fit(hyperbola,
                                                concentrations,
                                                signals,
                                                bounds=((yint_min, delta_y_min, sys.float_info.min*10000),
                                                        (yint_max, delta_y_max, np.inf)),
                                                sigma=sigmas)
    yint_stddev = covariance[0, 0] ** 0.5
    delta_y_stddev = covariance[1, 1] ** 0.5
    kd_stddev = covariance[2, 2] ** 0.5
    return yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev


def fit_quadratic(concentrations, signals):
    """
    Fit a titration to a quadratic.

    :param concentrations: X-axis values representing concentrations in arbitrary units
    :param signals: Y-axis values representing some kind of signal. Don't normalize this.

    :return:
        yint: the Y-intercept of the fit, often the background signal
        yint_stddev: standard deviation of the error of yint
        delta_y: total height of the fit
        delta_y_stddev: standard deviation of the error of delta_y
        kd: the dissociation constant
        kd_stddev: the standard deviation of the error of kd
        constant: the concentration of the substance that is not being titrated
        constant_stddev: the standard deviation of the error of constant

    """
    slope, intercept, _, _, _ = stats.linregress(concentrations, signals)
    yint_min, yint_max = (0.0, np.inf) if intercept > 0 else (-np.inf, 0.0)
    delta_y_min, delta_y_max = (0.0, np.inf) if slope > 0 else (-np.inf, 0.0)
    sigmas = [np.std(s) for s in signals]
    if 0 in sigmas:
        # if we have scalar values, the standard deviation will be zero everywhere and this will cause a division by
        # zero to occur. In this case, we just weight all points equally
        sigmas = np.ones((len(signals),))
    (yint, delta_y, kd, constant), covariance = curve_fit(quadratic,
                                                          concentrations,
                                                          signals,
                                                          bounds=((yint_min, delta_y_min, sys.float_info.min*10000, 0.0),
                                                                  (yint_max, delta_y_max, np.inf, np.inf)),
                                                          sigma=sigmas)
    yint_stddev = covariance[0, 0] ** 0.5
    delta_y_stddev = covariance[1, 1] ** 0.5
    kd_stddev = covariance[2, 2] ** 0.5
    constant_stddev = covariance[3, 3] ** 0.5
    return yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev, constant, constant_stddev
