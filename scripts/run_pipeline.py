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

from compostaje import bootstrap, data, fitting, pipeline, plots  # noqa: E402

RESULTS = ROOT / "results"
FIGS = ROOT / "docs" / "figuras"


def main(n_boot: int) -> None:
    RESULTS.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    raw = data.load_raw(data.download_raw())
    print(f"Filas leídas: {len(raw)}")

    res = pipeline.fit_dataset(raw)
    res.to_csv(RESULTS / "ajustes_por_serie.csv")
    counts = pipeline.winner_counts(res)
    counts.to_csv(RESULTS / "ganadores_por_criterio.csv")
    print(f"Series ajustadas: {len(res)}\n{counts}")

    # Bootstrap de las series representativas, con la numeración consistente.
    allser = data.assign_series_ids(raw)
    rows = []
    for uid in pipeline.BOOTSTRAP_IDS:
        s = data.extract_series(allser[allser["Unique_ID"] == uid], require_temp=False)
        b = bootstrap.bootstrap_series(s.t, s.c, n_boot=n_boot, seed=42)
        rows.append({"ID": uid, "Grupo_Residuo": s.group, "N": s.n,
                     **{f"{p}_{q}": b[p][q] for p in ("k", "Rm", "lambda")
                        for q in ("mean", "lo", "hi", "n_ok")}})
    boot = pd.DataFrame(rows)
    boot.to_csv(RESULTS / "bootstrap_series_representativas.csv", index=False)

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
