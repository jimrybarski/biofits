import numpy as np
from biofits import hyperbola, quadratic
import pytest


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

