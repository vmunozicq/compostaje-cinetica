"""Incertidumbre paramétrica por bootstrap (remuestreo de pares con reemplazo).

Reproduce ``analisis_bootstrap.py`` (V3) de la memoria: B = 1000 réplicas,
re-ajuste de primer orden (k) y Gompertz (Rm, λ), IC 95 % por percentiles.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from .models import first_order, gompertz_degraded

# Cotas usadas en el script de bootstrap original (más estrechas que en el ajuste
# principal; ver docs/supuestos.md).
BOOT_BOUNDS = {"k_max": 5.0, "rm_max": 50.0, "lam_max": 30.0}


def bootstrap_series(t, c, n_boot: int = 1000, seed: int | None = 42,
                     bounds: dict = BOOT_BOUNDS) -> dict:
    t = np.asarray(t, float)
    c = np.asarray(c, float)
    c0 = c[0]
    rng = np.random.default_rng(seed)
    ks, rms, lams = [], [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (OptimizeWarning, RuntimeWarning))
        for _ in range(n_boot):
            idx = rng.integers(0, len(t), len(t))
            tb, cb = t[idx], c[idx]
            yb = c0 - cb
            try:
                popt, _ = curve_fit(lambda tt, k, ci: first_order(tt, k, ci, c0), tb, cb,
                                    p0=[0.1, c[-1]], bounds=([0, 0], [bounds["k_max"], c0]),
                                    maxfev=1000)
                ks.append(popt[0])
            except Exception:  # noqa: BLE001 - réplica descartada, como en el original
                pass
            try:
                popt, _ = curve_fit(gompertz_degraded, tb, yb, p0=[(c0 - c).max(), 1, 1],
                                    bounds=([0, 0, 0], [c0, bounds["rm_max"], bounds["lam_max"]]),
                                    maxfev=1000)
                rms.append(popt[1])
                lams.append(popt[2])
            except Exception:  # noqa: BLE001
                pass

    def summary(x):
        x = np.asarray(x)
        if len(x) == 0:
            return dict(mean=np.nan, lo=np.nan, hi=np.nan, n_ok=0, samples=x)
        lo, hi = np.percentile(x, [2.5, 97.5])
        return dict(mean=float(x.mean()), lo=float(lo), hi=float(hi), n_ok=len(x), samples=x)

    return {"k": summary(ks), "Rm": summary(rms), "lambda": summary(lams)}
