import numpy as np
from biofits import hyperbola, quadratic, fit_hyperbola, fit_quadratic
import pytest
from hypothesis import given, strategies as st
np.seterr(all='raise')

# Some test values to play with
concentrations1 = np.array([0, 10, 50, 100, 150, 200, 250, 280, 400])
fluorescence1 = np.array([0.0, 0.14, 0.30, 0.45, 0.66, 0.75, 0.93, 0.99, 1.0])


def test_hyperbola():
    concentrations = np.array([0, 10, 90])
    f0 = 0.0
    delta_f = 1.0
    kd = 10.0
    first, second, third = hyperbola(concentrations, f0, delta_f, kd)
    assert first == 0.0
    assert second == 0.5
    assert third == 0.9


def test_hyperbola_f0():
    concentrations = np.array([0, 10, 90])
    delta_f = 1.0
    kd = 10.0
    a, b, c = hyperbola(concentrations, 0.0, delta_f, kd)
    d, e, f = hyperbola(concentrations, 3.8, delta_f, kd)
    for first, second in ((a, d), (b, e), (c, f)):
        assert second - first == pytest.approx(3.8)


def test_delta_f():
    concentrations = np.array([1.0, 10.0, 90.0])
    f0 = 0.0
    kd = 0.0
    for delta_f in (-30, -2.0, -1.0, 0.0, 10.0, 50.0, 1000.0, 10000000.0):
        a, b, c = hyperbola(concentrations, f0, delta_f, kd)
        assert a == delta_f
        assert b == delta_f
        assert c == delta_f


def test_kd():
    f0 = 0.0
    delta_f = 1.0
    kd = 100.0
    for i in np.linspace(0.0, 10000000000000.0, num=1000):
        concentrations = np.array([0.0, 10.0, i])
        a, b, c = hyperbola(concentrations, f0, delta_f, kd)
        assert c < 1


def test_limit():
    f0 = 0.0
    delta_f = 1.0
    kd = 0.0000001
    for i in np.linspace(20.0, 10000000000000.0, num=1000):
        concentrations = np.array([0.0, 10.0, i])
        a, b, c = hyperbola(concentrations, f0, delta_f, kd)
        assert c < i
        assert c == pytest.approx(1.0)


def test_hyperbolic_fit():
    f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev = fit_hyperbola(concentrations1, fluorescence1)
    assert kd == pytest.approx(251.01, rel=0.001)
    assert kd_stddev == pytest.approx(84.44, rel=0.001)


def test_hyperbolic_fit_inhibition():
    f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev = fit_hyperbola(concentrations1, fluorescence1[::-1])
    assert kd > 0
    assert kd_stddev > 0
    assert delta_f < 0
    assert f0 > 0


def test_quadratic_fit():
    f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev, constant, constant_stddev = fit_quadratic(concentrations1, fluorescence1)
    assert kd < 200
    assert constant > 1


def test_quadratic_fit_inhibition():
    f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev, constant, constant_stddev = fit_quadratic(concentrations1, fluorescence1[::-1])
    assert kd > 0
    assert kd_stddev > 0
    assert delta_f < 0
    assert f0 > 0


@given(st.lists(st.floats(allow_nan=True, allow_infinity=True), min_size=3), st.lists(st.floats(allow_nan=True, allow_infinity=True), min_size=3))
def test_hyperbolic_fit_hypothesis(concentrations, fluorescence):
    concentrations = np.array(sorted(concentrations))
    fluorescence = np.array(sorted(fluorescence))
    if len(concentrations) != len(fluorescence) or len(set(concentrations)) != len(fluorescence) or len(concentrations) != len(set(fluorescence)) or len(concentrations) < 3 or len(fluorescence) < 3 or np.any(np.isnan(fluorescence)) or np.any(np.isnan(concentrations)) or len(set(concentrations)) != len(concentrations) or not any(fluorescence) or np.any(concentrations[concentrations < 10**-50]) or np.any(concentrations[concentrations > 10**50]) or np.max(fluorescence) - np.min(fluorescence) < 10**-50:
        with pytest.raises(ValueError):
            fit_hyperbola(concentrations, fluorescence)
    else:
        f0, f0_stddev, delta_f, delta_f_stddev, kd, kd_stddev = fit_hyperbola(concentrations, fluorescence)
        assert kd > 0.0
