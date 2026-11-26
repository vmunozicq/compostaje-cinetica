import numpy as np
import pytest

from compostaje import fitting, models

T = np.array([0, 2, 4, 7, 10, 14, 21, 28, 35], float)


def test_recupera_parametros_de_primer_orden():
    c = models.first_order(T, 0.15, 20, 50)
    r = fitting.fit_first_order(T, c)
    assert r.params["k"] == pytest.approx(0.15, rel=1e-4)
    assert r.params["c_inf"] == pytest.approx(20, rel=1e-4)
    assert r.r2 == pytest.approx(1)


def test_recupera_gompertz_con_latencia():
    c = models.gompertz(T, 25, 3, 5, 50)
    r = fitting.fit_gompertz(T, c)
    assert r.params["lambda"] == pytest.approx(5, rel=1e-3)
    assert r.params["Rm"] == pytest.approx(3, rel=1e-3)


def test_aicc_prefiere_el_modelo_simple_si_basta():
    rng = np.random.default_rng(0)
    c = models.first_order(T, 0.2, 15, 45) + rng.normal(0, 0.3, len(T))
    fits = fitting.fit_all(T, c)
    assert fitting.best_model(fits) == "first_order"


def test_bifasico_etiqueta_k1_como_rapida():
    c = models.biphasic(T, 0.05, 0.6, 0.4, 10, 50)
    r = fitting.fit_biphasic(T, c)
    assert r.params["k1"] > r.params["k2"]


def test_modo_consistente_mide_el_tiempo_desde_el_primer_dato():
    """Serie que parte en t = 5 d: solo el modo consistente recupera k (L2)."""
    t = T + 5
    c = models.first_order(T, 0.15, 20, 50)  # la cinética empieza en el primer dato
    r = fitting.fit_all(t, c, modo="consistente")["first_order"]
    assert r.params["k"] == pytest.approx(0.15, rel=1e-4)
    assert r.predict(t) == pytest.approx(c, abs=1e-6)  # predict usa el tiempo original
    r_mem = fitting.fit_all(t, c, modo="memoria")["first_order"]
    assert r_mem.r2 < r.r2


def test_modo_invalido():
    with pytest.raises(ValueError):
        fitting.fit_all(T, models.first_order(T, 0.1, 20, 50), modo="otro")


def test_arrhenius_tref_recupera_parametros():
    temp = np.array([25, 45, 60, 62, 55, 48, 40, 33, 28], float)
    t_ref = temp.mean()
    c = models.arrhenius_tref(T, 0.12, 60000, 25, 50, T, temp, t_ref)
    r = fitting.fit_arrhenius_tref(T, c, temp)
    assert r.params["T_ref"] == pytest.approx(t_ref)
    assert r.params["k_ref"] == pytest.approx(0.12, rel=1e-2)
    assert r.params["Ea"] == pytest.approx(60000, rel=2e-2)
    assert r.r2 == pytest.approx(1, abs=1e-6)
