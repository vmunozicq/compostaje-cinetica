"""Reproduce todos los resultados y figuras del repositorio.

Uso:  python scripts/run_pipeline.py [--boot 1000]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from compostaje import bootstrap, data, fitting, ml, pipeline, plots, reactor  # noqa: E402

RESULTS = ROOT / "results"
FIGS = ROOT / "docs" / "figuras"


def main(n_boot: int) -> None:
    RESULTS.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    raw = data.load_raw(data.download_raw())
    print(f"Filas leídas: {len(raw)}")

    # Versión corregida (t − t0, Arrhenius con T_ref) y reproducción de la memoria.
    res = pipeline.fit_dataset(raw, modo="consistente")
    res.to_csv(RESULTS / "ajustes_por_serie.csv")
    res_mem = pipeline.fit_dataset(raw, modo="memoria")
    res_mem.to_csv(RESULTS / "ajustes_modo_memoria.csv")
    counts = pipeline.winner_counts(res, res_mem)
    counts.to_csv(RESULTS / "ganadores_por_criterio.csv")
    print(f"Series ajustadas: {len(res)}\n{counts}")

    # Bootstrap de las series representativas, con la numeración consistente.
    allser = data.assign_series_ids(raw)
    boots = {}
    for nombre, metodo, modo in (("bootstrap_series_representativas", "residuos", "consistente"),
                                 ("bootstrap_modo_memoria", "pares", "memoria")):
        rows = []
        for uid in pipeline.BOOTSTRAP_IDS:
            s = data.extract_series(allser[allser["Unique_ID"] == uid], require_temp=False)
            b = bootstrap.bootstrap_series(s.t, s.c, n_boot=n_boot, seed=42,
                                           metodo=metodo, modo=modo)
            rows.append({"ID": uid, "Grupo_Residuo": s.group, "N": s.n,
                         **{f"{p}_{q}": b[p][q] for p in ("k", "Rm", "lambda")
                            for q in ("mean", "lo", "hi", "n_ok")}})
        boots[nombre] = pd.DataFrame(rows)
        boots[nombre].to_csv(RESULTS / f"{nombre}.csv", index=False)
    boot = boots["bootstrap_series_representativas"]

    # Aprendizaje automático: ¿se puede predecir k sin ensayo? (docs/modelo_ml.md)
    tab = ml.build_ml_table(raw, res)
    tab.to_csv(RESULTS / "ml_tabla.csv")
    cv = ml.cross_validate(tab)
    cv.to_csv(RESULTS / "ml_validacion.csv", index=False)
    pred = ml.evaluar_dominio(tab)
    pred.to_csv(RESULTS / "ml_predicciones_cv.csv")
    print(cv.round(3).to_string(index=False))
    plots.plot_ml_paridad(pred, FIGS / "ml_paridad.png")

    # Simulador de reactor (docs/modelo_reactor.md)
    exp = reactor.cargar_experimento(ROOT / "data" / "reference"
                                     / "petric2008_exp2_digitalizado.csv")
    sims = {"Simulador (parámetros del artículo)": reactor.simular(),
            "Con todas las correcciones de la Ec. 2": reactor.simular(
                par=reactor.Parametros(correcciones="todas"))}
    val = pd.concat([reactor.validar(s, exp).assign(variante=k) for k, s in sims.items()])
    val.to_csv(RESULTS / "reactor_validacion.csv", index=False)
    print(val.round(2).to_string(index=False))
    plots.plot_reactor_validacion(sims, exp, FIGS / "reactor_validacion.png")
    reactor.variar("caudal_m3_h", [0.06, 0.09, 0.10, 0.14, 0.18, 0.22, 0.26]).to_csv(
        RESULTS / "reactor_caudal.csv", index=False)

    # Figuras
    plots.plot_winners(counts, FIGS / "ganadores_aicc.png")
    plots.plot_r2_box(res, FIGS / "r2_por_modelo.png")
    plots.plot_k_by_group(res, FIGS / "k_por_sustrato.png")
    plots.plot_bootstrap(boot, FIGS / "bootstrap_k.png")
    plots.plot_k_vs(res, "Init_MC", "Humedad inicial (%)", FIGS / "k_vs_humedad.png")
    plots.plot_k_vs(res, "Init_CN", "C/N inicial", FIGS / "k_vs_cn.png")

    valid = data.build_validated_dataset(raw)
    for uid, label in pipeline.CASE_STUDY_IDS.items():
        s = data.extract_series(valid[valid["Unique_ID"] == uid])
        fits = fitting.fit_all(s.t, s.c, s.temp)
        best = pipeline.SHORT[fitting.best_model(fits)]
        fig, ax = plt.subplots(figsize=(7, 4.2))
        plots.plot_series_fits(s, fits, ax, title=f"Serie {uid} · {label} · mejor por AICc: "
                                                  f"{plots.SHORT_LABELS[best]}")
        fig.tight_layout()
        fig.savefig(FIGS / f"caso_{uid}.png", dpi=200)
        plt.close(fig)
    print("Listo: results/ y docs/figuras/ actualizados.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000, help="réplicas bootstrap por serie")
    main(ap.parse_args().boot)
