"""Tests del dimensionamiento de planta (docs/modelo_diseno.md)."""

import math
from pathlib import Path

import pandas as pd
import pytest

from compostaje import diseno

RES = Path(__file__).resolve().parents[1] / "results" / "ajustes_por_serie.csv"


@pytest.fixture(scope="module")
def raz():
    return diseno.razones_reactor()


@pytest.fixture(scope="module")
def cin():
    return diseno.cinetica_sustrato(pd.read_csv(RES), "Kitchen Waste")


def test_tiempo_residencia_primer_orden():
    # Para X = 1 − e^(−1), t = 1/k.
    assert diseno.tiempo_residencia_d(0.1, 1 - math.exp(-1)) == pytest.approx(10.0)
    with pytest.raises(ValueError):
        diseno.tiempo_residencia_d(0.1, 1.0)
    with pytest.raises(ValueError):
        diseno.tiempo_residencia_d(0.0, 0.5)


def test_cinetica_sustrato_ordenada(cin):
    assert cin.n_series == 45  # series de cocina con R² ≥ 0,5
    assert cin.k_p25_d < cin.k_p50_d < cin.k_p75_d
    assert 0 < cin.frac_degradable < 1


def test_razones_del_reactor_validado(raz):
    # Aire del experimento 2: 0,18 m³/h × 336 h sobre ~1,26 kg de MO degradada.
    assert raz.aire_m3_por_kg_om == pytest.approx(47.9, rel=0.02)
    assert raz.evaporacion_kg_por_kg_om == pytest.approx(2.75, rel=0.02)


def test_dimensionar_escala_y_conservador(cin, raz):
    ent = diseno.EntradaPlanta(10, 0.6, 0.8, 0.9, 600, 2.0)
    t = diseno.dimensionar(ent, cin, raz)
    assert list(t["escenario"]) == ["conservador (P25)", "típico (P50)", "optimista (P75)"]
    assert t["volumen_m3"].is_monotonic_decreasing  # k bajo → más tiempo y volumen
    doble = diseno.dimensionar(diseno.EntradaPlanta(20, 0.6, 0.8, 0.9, 600, 2.0), cin, raz)
    assert (doble["volumen_m3"] / t["volumen_m3"]).round(9).eq(2).all()
    aire = "aire_control_termico_m3_h"
    assert (doble[aire] / t[aire]).round(9).eq(2).all()
    # Balance: O₂ = y_O₂ · MO degradada.
    assert t["o2_kg_d"].iloc[0] == pytest.approx(0.705 * t["om_degradada_kg_d"].iloc[0])


def test_k_medido_por_el_usuario(cin, raz):
    ent = diseno.EntradaPlanta(10, 0.6, 0.8, 0.9, 600, 2.0)
    t = diseno.dimensionar(ent, cin, raz, k_usuario_d=0.1)
    assert len(t) == 1
    assert t["tiempo_residencia_d"].iloc[0] == pytest.approx(-math.log(0.1) / 0.1)


def test_entradas_invalidas():
    with pytest.raises(ValueError):
        diseno.EntradaPlanta(-1, 0.6, 0.8, 0.9, 600, 2.0)
    with pytest.raises(ValueError):
        diseno.EntradaPlanta(10, 1.2, 0.8, 0.9, 600, 2.0)
    with pytest.raises(ValueError):
        diseno.EntradaPlanta(10, 0.6, 0.8, 0.9, 0, 2.0)
