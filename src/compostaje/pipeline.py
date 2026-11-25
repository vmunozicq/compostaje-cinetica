"""Flujo completo: datos → series → ajuste de 4 modelos → selección por AICc."""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import data, fitting, metrics
from .models import N_PARAMS

MODELS = ["first_order", "biphasic", "gompertz", "arrhenius"]
SHORT = {"first_order": "1st", "biphasic": "Bi", "gompertz": "Gomp", "arrhenius": "Arr"}

# Series usadas en la memoria para el bootstrap (Figura 10) y los casos de estudio.
BOOTSTRAP_IDS = ["8_1", "27_1", "3_1", "25_1", "461_1", "459_1", "152_1", "161_1", "454_1", "404_1"]
CASE_STUDY_IDS = {"25_1": "Residuos de cocina", "461_1": "Lodos sanitarios", "152_1": "Estiércol"}


def raw_series_stats(df_raw: pd.DataFrame) -> pd.DataFrame:
    """N y SST de cada serie *antes* del filtro físico-biológico.

    La memoria calculó AIC/AICc/BIC con estos N y SST (ver docs/supuestos.md, D1).
    """
    allser = data.assign_series_ids(df_raw)
    rows = {}
    for uid, g in allser.groupby("Unique_ID"):
        var = "TOC_percent" if g["TOC_percent"].notna().sum() >= data.MIN_POINTS else "OM_percent"
        y = g[["Time", var]].dropna()[var].to_numpy(float)
        if len(y):
            rows[uid] = {"N_raw": len(y), "SST_raw": float(np.sum((y - y.mean()) ** 2))}
    return pd.DataFrame.from_dict(rows, orient="index")


def fit_dataset(df_raw: pd.DataFrame) -> pd.DataFrame:
    """Ajusta los 4 modelos a todas las series válidas y devuelve una fila por serie."""
    valid = data.build_validated_dataset(df_raw)
    rows = []
    for s in data.iter_series(valid):
        fits = fitting.fit_all(s.t, s.c, s.temp)
        row = {"ID": s.uid, "Grupo_Residuo": s.group, "Residuo_Original": s.waste_type,
               "Variable": s.variable, "N": s.n, "SST": metrics.sst(s.c),
               "Init_MC": s.init_mc, "Init_CN": s.init_cn}
        for r in fits.values():
            row.update(r.as_row())
        for crit in ("aic", "aicc", "bic"):
            best = fitting.best_model(fits, crit)
            row[f"best_{crit}"] = SHORT.get(best, "Ninguno")
        rows.append(row)
    out = pd.DataFrame(rows).set_index("ID")
    return add_legacy_criteria(out, raw_series_stats(df_raw))


def add_legacy_criteria(res: pd.DataFrame, raw_stats: pd.DataFrame) -> pd.DataFrame:
    """Recalcula AICc como en la memoria: SSE = (1 − R²)·SST_raw y n = N_raw."""
    res = res.join(raw_stats)
    for m in MODELS:
        vals = []
        for _, row in res.iterrows():
            r2 = row.get(f"{m}__r2")
            if pd.isna(r2) or r2 <= 0 or row["SST_raw"] <= 0:
                vals.append(np.inf)
                continue
            vals.append(metrics.information_criteria((1 - r2) * row["SST_raw"],
                                                     int(row["N_raw"]), N_PARAMS[m])[1])
        res[f"{m}__aicc_legacy"] = vals
    cols = [f"{m}__aicc_legacy" for m in MODELS]
    def pick(row):
        finite = row[np.isfinite(row.astype(float))]
        return SHORT[finite.astype(float).idxmin().split("__")[0]] if len(finite) else "Ninguno"

    res["best_aicc_legacy"] = res[cols].apply(pick, axis=1)
    return res


def winner_counts(res: pd.DataFrame) -> pd.DataFrame:
    order = ["1st", "Gomp", "Bi", "Arr", "Ninguno"]
    out = pd.DataFrame({
        "AICc (n y SST consistentes)": res["best_aicc"].value_counts(),
        "AICc (cálculo de la memoria)": res["best_aicc_legacy"].value_counts(),
        "BIC": res["best_bic"].value_counts(),
    }).reindex(order).fillna(0).astype(int)
    return out.loc[out.sum(axis=1) > 0]
