from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from compostaje import reactor

REF = Path(__file__).resolve().parents[1] / "data" / "reference"
EXP = REF / "petric2008_exp2_digitalizado.csv"


@pytest.fixture(scope="module")
def sim():
    return reactor.simular()


def test_constantes_termodinamicas():
    assert reactor.p_sat(373.15) == pytest.approx(101325, rel=0.02)  # ebullición
    assert reactor.henry("O2", 298.15) == pytest.approx(4.4e9, rel=0.1)
    assert reactor.henry("CO2", 298.15) == pytest.approx(1.65e8, rel=0.1)


def test_balance_de_masa_cierra(sim):
    """La materia orgánica solo baja y el O₂ de salida nunca supera al del aire."""
    assert (np.diff(sim["m_OM_kg"]) <= 1e-9).all()
    assert sim["O2_pct"].max() <= 21.0 + 1e-6


def test_valida_contra_experimento_2(sim):
    """Criterios (docs/modelo_reactor.md): conversión RMSE < 3 puntos, T RMSE < 7 °C,
    O₂ RMSE < 3,5 % vol. El modelo del propio artículo tiene error medio de 2,3 °C y
    error máximo de 14,3 °C en temperatura (Tabla 6)."""
    v = reactor.validar(sim, reactor.cargar_experimento(EXP)).set_index("variable")
    assert v.loc["conversion_pct", "RMSE"] < 3
    assert v.loc["T_C", "RMSE"] < 7
    assert v.loc["O2_pct", "RMSE"] < 3.5


def test_mas_aire_enfria():
    v = reactor.variar("caudal_m3_h", [0.09, 0.26])
    assert v["T_max_C"].iloc[0] > v["T_max_C"].iloc[1]


def test_entradas_invalidas():
    with pytest.raises(ValueError):
        reactor.Operacion(humedad=1.2)
    with pytest.raises(ValueError):
        reactor.simular(par=replace(reactor.Parametros(), correcciones="otra"))
