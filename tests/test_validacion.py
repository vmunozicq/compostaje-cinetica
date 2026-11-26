"""Validación contra los resultados originales de la memoria de título.

Referencias en data/reference/ (salidas de los scripts originales de 2025):
- memoria_tabla_criterios.csv: R², N, SST y modelo ganador de las 159 series.

Se ajusta en modo "memoria" (sin desplazar t0 y con Arrhenius en la forma original),
que es el que reproduce los scripts de 2025.

Tolerancias:
- R² de primer orden, bifásico y Gompertz: 1e-4. El algoritmo (curve_fit, mismas
  cotas y valores iniciales) es idéntico; solo cambian versiones de SciPy/NumPy.
- Arrhenius no se compara punto a punto: L-BFGS-B con gradiente numérico sobre un
  ODE converge a mínimos locales distintos según la versión (ver docs/supuestos.md).
"""

from pathlib import Path

import pandas as pd
import pytest

from compostaje import data, pipeline

REF = Path(__file__).resolve().parents[1] / "data" / "reference" / "memoria_tabla_criterios.csv"
TOL_R2 = 1e-4


@pytest.fixture(scope="module")
def resultados():
    raw = data.load_raw(data.DEFAULT_RAW_PATH)
    return pipeline.fit_dataset(raw, modo="memoria"), pd.read_csv(REF).set_index("ID")


def test_mismas_159_series(resultados):
    res, ref = resultados
    assert len(res) == 159
    assert set(res.index) == set(ref.index)


@pytest.mark.parametrize("modelo,col_ref", [("first_order", "R2_1st"),
                                            ("biphasic", "R2_Bi"),
                                            ("gompertz", "R2_Gomp")])
def test_r2_identico_a_la_memoria(resultados, modelo, col_ref):
    res, ref = resultados
    r2_ref = ref.loc[res.index, col_ref]
    # -999 marca un ajuste que no convergió en el script original. Si una versión más
    # nueva de NumPy/SciPy sí converge, no es una discrepancia del modelo: se excluye.
    ok = r2_ref > -999
    diff = (res.loc[ok, f"{modelo}__r2"] - r2_ref[ok]).abs()
    assert diff.max() < TOL_R2


def test_k_primer_orden_identica(resultados):
    res, ref = resultados
    diff = (res["first_order__k"] - ref.loc[res.index, "k_1st"]).abs()
    assert diff.max() < 1e-5


def test_modo_memoria_reproduce_ganadores_aicc(resultados):
    """Con N y SST sin filtrar (cálculo original) se recuperan los ganadores de la memoria.

    Única diferencia: 2 series sin ningún ajuste con R² > 0, que el script original
    rotulaba "1st" por defecto y aquí quedan como "Ninguno".
    """
    res, ref = resultados
    distintos = res.index[res["best_aicc_legacy"] != ref.loc[res.index, "Best_AICc"]]
    assert set(distintos) <= {"4_1", "188_1"}
    assert (res.loc[distintos, "best_aicc_legacy"] == "Ninguno").all()
