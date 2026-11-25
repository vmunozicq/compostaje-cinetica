import numpy as np
import pytest

from compostaje import metrics


def test_r2_perfecto_y_media():
    y = np.array([1.0, 2.0, 4.0, 8.0])
    assert metrics.r2(y, y) == pytest.approx(1)
    assert metrics.r2(y, np.full_like(y, y.mean())) == pytest.approx(0)


def test_criterios_de_informacion_a_mano():
    sse, n, p = 12.0, 8, 2
    aic, aicc, bic = metrics.information_criteria(sse, n, p)
    base = 8 * np.log(12 / 8)
    assert aic == pytest.approx(base + 4)
    assert aicc == pytest.approx(base + 4 + 2 * 2 * 3 / 5)
    assert bic == pytest.approx(base + 2 * np.log(8))


def test_aicc_indefinido_con_pocos_datos():
    assert metrics.information_criteria(1.0, 5, 4) == (np.inf, np.inf, np.inf)


def test_nrmse_en_porcentaje():
    y = np.array([10.0, 10.0])
    assert metrics.nrmse_pct(y, y + 1) == pytest.approx(10)
