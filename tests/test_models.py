import numpy as np
import pytest

from compostaje import models


def test_primer_orden_condiciones_borde():
    assert models.first_order(0, 0.2, 10, 50) == pytest.approx(50)
    assert models.first_order(1e4, 0.2, 10, 50) == pytest.approx(10)


def test_bifasico_se_reduce_a_primer_orden():
    t = np.linspace(0, 30, 20)
    np.testing.assert_allclose(models.biphasic(t, 0.1, 0.1, 0.3, 5, 40),
                               models.first_order(t, 0.1, 5, 40))


def test_gompertz_latencia_y_asintota():
    # Antes de λ casi no hay degradación; a tiempo largo se degrada P.
    assert models.gompertz_degraded(0, 20, 2, 10) < 0.01
    assert models.gompertz_degraded(1e3, 20, 2, 10) == pytest.approx(20)


def test_gompertz_tasa_maxima_es_rm():
    p, rm, lam = 30.0, 2.5, 4.0
    t = np.linspace(0, 60, 60001)
    pendiente = np.gradient(models.gompertz_degraded(t, p, rm, lam), t)
    assert pendiente.max() == pytest.approx(rm, rel=1e-3)


def test_arrhenius_isotermico_igual_a_primer_orden():
    t = np.linspace(0, 20, 11)
    temp = np.full_like(t, 55.0)
    a, ea = 5e4, 40000.0
    k = models.k_arrhenius(55.0, a, ea)
    np.testing.assert_allclose(models.arrhenius(t, a, ea, 8, 45, t, temp),
                               models.first_order(t, k, 8, 45), rtol=1e-5)
