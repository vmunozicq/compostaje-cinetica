import numpy as np
import pandas as pd
import pytest

from compostaje import ml


@pytest.fixture
def tabla_sintetica():
    """k depende solo de la temperatura media: un buen modelo debe superar a la mediana."""
    rng = np.random.default_rng(0)
    n = 90
    tab = pd.DataFrame({
        "Grupo_Residuo": rng.choice(["Kitchen Waste", "Manure", "Sludge"], n),
        "Init_pH": rng.uniform(6, 8, n), "Init_CN": rng.uniform(15, 35, n),
        "Init_MC": rng.uniform(50, 70, n), "T_media": rng.uniform(25, 55, n),
    }, index=[f"{i}_1" for i in range(n)])
    tab[ml.TARGET] = -2 + 0.03 * tab["T_media"] + rng.normal(0, 0.05, n)
    tab["grupo_cv"] = np.arange(n)
    return tab


def test_modelo_supera_a_la_mediana_cuando_hay_senal(tabla_sintetica):
    cv = ml.cross_validate(tabla_sintetica).set_index(["Validación", "Modelo"])
    base = cv.loc[("Series nuevas", "Mediana (referencia)"), "RMSE_log10"]
    assert cv.loc[("Series nuevas", "Ridge"), "RMSE_log10"] < 0.5 * base
    assert cv.loc[("Series nuevas", "Ridge"), "R2"] > 0.9


def test_dominio_marca_entradas_lejanas(tabla_sintetica):
    dom = ml.DominioAplicabilidad.fit(tabla_sintetica)
    centro = tabla_sintetica[ml.NUM_FEATURES].median().to_dict()
    assert dom.evaluar(centro)["dentro"]
    lejos = dict(centro, Init_CN=80)
    ev = dom.evaluar(lejos)
    assert not ev["dentro"] and ev["fuera_de_rango"] == ["Init_CN"]


def test_predict_k_valida_humedad(tabla_sintetica):
    model, _ = ml.fit_final(tabla_sintetica, "Ridge")
    entrada = dict(tabla_sintetica[ml.NUM_FEATURES].median().to_dict(), Grupo_Residuo="Manure")
    assert ml.predict_k(model, entrada) > 0
    with pytest.raises(ValueError):
        ml.predict_k(model, dict(entrada, Init_MC=120))
