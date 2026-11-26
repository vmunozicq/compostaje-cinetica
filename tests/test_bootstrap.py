import numpy as np
import pytest

from compostaje import bootstrap, models

T = np.array([0, 2, 4, 7, 10, 14, 21, 28], float)


@pytest.mark.parametrize("metodo", ["residuos", "pares"])
def test_ic_contiene_el_k_verdadero(metodo):
    rng = np.random.default_rng(1)
    c = models.first_order(T, 0.2, 20, 50) + rng.normal(0, 0.4, len(T))
    b = bootstrap.bootstrap_series(T, c, n_boot=200, seed=0, metodo=metodo)
    assert b["k"]["lo"] < 0.2 < b["k"]["hi"]
    assert b["k"]["n_ok"] > 150


def test_misma_semilla_mismo_resultado():
    c = models.first_order(T, 0.2, 20, 50) + np.linspace(-0.3, 0.3, len(T))
    b1 = bootstrap.bootstrap_series(T, c, n_boot=50, seed=7)
    b2 = bootstrap.bootstrap_series(T, c, n_boot=50, seed=7)
    assert b1["k"]["lo"] == b2["k"]["lo"]


def test_metodo_invalido():
    with pytest.raises(ValueError):
        bootstrap.bootstrap_series(T, T, metodo="jackknife")
