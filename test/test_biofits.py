import numpy as np
from biofits import hyperbola, quadratic, fit_hyperbola, fit_quadratic
import pytest

# Some test values to play with
concentrations1 = np.array([0, 10, 50, 100, 150, 200, 250, 280, 400])
fluorescence1 = np.array([0.0, 0.14, 0.30, 0.45, 0.66, 0.75, 0.93, 0.99, 1.0])


def test_hyperbola():
    concentrations = np.array([0, 10, 90])
    yint = 0.0
    delta_y = 1.0
    kd = 10.0
    first, second, third = hyperbola(concentrations, yint, delta_y, kd)
    assert first == 0.0
    assert second == 0.5
    assert third == 0.9


def test_hyperbola_yint():
    concentrations = np.array([0, 10, 90])
    delta_y = 1.0
    kd = 10.0
    a, b, c = hyperbola(concentrations, 0.0, delta_y, kd)
    d, e, f = hyperbola(concentrations, 3.8, delta_y, kd)
    for first, second in ((a, d), (b, e), (c, f)):
        assert second - first == pytest.approx(3.8)


def test_delta_y():
    concentrations = np.array([1.0, 10.0, 90.0])
    yint = 0.0
    kd = 0.0
    for delta_y in (-30, -2.0, -1.0, 0.0, 10.0, 50.0, 1000.0, 10000000.0):
        a, b, c = hyperbola(concentrations, yint, delta_y, kd)
        assert a == delta_y
        assert b == delta_y
        assert c == delta_y


def test_kd():
    yint = 0.0
    delta_y = 1.0
    kd = 100.0
    for i in np.linspace(0.0, 10000000000000.0, num=1000):
        concentrations = np.array([0.0, 10.0, i])
        a, b, c = hyperbola(concentrations, yint, delta_y, kd)
        assert c < 1


def test_limit():
    yint = 0.0
    delta_y = 1.0
    kd = 0.0000001
    for i in np.linspace(20.0, 10000000000000.0, num=1000):
        concentrations = np.array([0.0, 10.0, i])
        a, b, c = hyperbola(concentrations, yint, delta_y, kd)
        assert c < i
        assert c == pytest.approx(1.0)


def test_hyperbolic_fit():
    yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev = fit_hyperbola(concentrations1, fluorescence1)
    assert kd == pytest.approx(251.01, rel=0.02)
    assert kd_stddev == pytest.approx(84.44, rel=0.02)


def test_hyperbolic_fit_inhibition():
    yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev = fit_hyperbola(concentrations1, fluorescence1[::-1])
    assert kd > 0
    assert kd_stddev > 0
    assert delta_y < 0
    assert yint > 0


def test_quadratic():
    concentrations = np.array([0, 10, 90000000])
    yint = 3.0
    delta_y = 2.0
    kd = 10.0
    first, second, third = quadratic(concentrations, yint, delta_y, kd, 1.0)
    assert first == yint
    assert (second - yint) / delta_y == pytest.approx(0.487507803)
    assert third - yint < delta_y


def test_quadratic_fit():
    yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev, constant, constant_stddev = fit_quadratic(concentrations1, fluorescence1)
    assert kd < 200
    assert constant > 1


def test_quadratic_fit_inhibition():
    yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev, constant, constant_stddev = fit_quadratic(concentrations1, fluorescence1[::-1])
    assert kd > 0
    assert kd_stddev > 0
    assert delta_y < 0
    assert yint > 0


def test_steves_real_data():
    concentrations = [0.001, 0.1, 0.3, 1, 3, 10, 30, 100, 300]
    signals = [0.297, 0.242, 0.353, 0.461, 0.543, 0.653, 0.780, 0.763, 0.701]
    yint, yint_stddev, delta_y, delta_y_stddev, kd, kd_stddev = fit_hyperbola(concentrations, signals)
    assert kd == pytest.approx(1.89, rel=0.05)
    assert kd_stddev == pytest.approx(0.66, rel=0.05)
