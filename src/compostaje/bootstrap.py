"""Incertidumbre paramétrica por bootstrap: k de primer orden, Rm y λ de Gompertz.

Dos métodos (docs/supuestos.md, L3):

- ``"pares"``: remuestreo de pares (t, C) con reemplazo. Es el método de la memoria
  (``analisis_bootstrap.py`` V3). Con 5 a 9 puntos produce réplicas con tiempos
  repetidos y ajustes degenerados.
- ``"residuos"``: se ajusta el modelo una vez, se remuestrean sus residuos y se suman a
  la curva ajustada, manteniendo los tiempos medidos. Es más estable con n pequeño,
  pero supone que el modelo es correcto y que los errores son homocedásticos
  (supuesto 5 de docs/modelo.md).

En ambos casos: re-ajuste de primer orden (k) y Gompertz (Rm, λ), IC 95 % por
percentiles. ``modo`` se interpreta igual que en :func:`fitting.fit_all` (t − t0).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from .models import first_order, gompertz_degraded

# Cotas usadas en el script de bootstrap original (más estrechas que en el ajuste
# principal; ver docs/supuestos.md).
BOOT_BOUNDS = {"k_max": 5.0, "rm_max": 50.0, "lam_max": 30.0}
METODOS = ("residuos", "pares")


def _fit_k(t, c, c0, bounds, p0_ci):
    popt, _ = curve_fit(lambda tt, k, ci: first_order(tt, k, ci, c0), t, c,
                        p0=[0.1, p0_ci], bounds=([0, 0], [bounds["k_max"], c0]), maxfev=1000)
    return popt


def _fit_gomp(t, y, c0, bounds, p0_p):
    popt, _ = curve_fit(gompertz_degraded, t, y, p0=[p0_p, 1, 1],
                        bounds=([0, 0, 0], [c0, bounds["rm_max"], bounds["lam_max"]]),
                        maxfev=1000)
    return popt


def bootstrap_series(t, c, n_boot: int = 1000, seed: int | None = 42,
                     bounds: dict = BOOT_BOUNDS, metodo: str = "residuos",
                     modo: str = "consistente") -> dict:
    """IC 95 % bootstrap de k (d⁻¹), Rm (%·d⁻¹) y λ (d) para una serie (t en d, C en %)."""
    if metodo not in METODOS:
        raise ValueError(f"metodo debe ser uno de {METODOS}, no {metodo!r}")
    if modo not in ("consistente", "memoria"):
        raise ValueError(f"modo debe ser 'consistente' o 'memoria', no {modo!r}")
    t = np.asarray(t, float)
    c = np.asarray(c, float)
    if modo == "consistente":
        t = t - t[0]
    c0 = c[0]
    y = c0 - c
    rng = np.random.default_rng(seed)
    ks, rms, lams = [], [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (OptimizeWarning, RuntimeWarning))
        if metodo == "residuos":
            try:
                k, ci = _fit_k(t, c, c0, bounds, c[-1])
                c_hat = first_order(t, k, ci, c0)
                res_c = c - c_hat
            except Exception:  # noqa: BLE001 - sin ajuste base no hay residuos
                c_hat = None
            try:
                p, rm, lam = _fit_gomp(t, y, c0, bounds, y.max())
                y_hat = gompertz_degraded(t, p, rm, lam)
                res_y = y - y_hat
            except Exception:  # noqa: BLE001
                y_hat = None
        for _ in range(n_boot):
            if metodo == "pares":
                idx = rng.integers(0, len(t), len(t))
                tb, cb, yb = t[idx], c[idx], y[idx]
            else:
                tb = t
                cb = None if c_hat is None else c_hat + rng.choice(res_c, len(t))
                yb = None if y_hat is None else y_hat + rng.choice(res_y, len(t))
            if cb is not None:
                try:
                    ks.append(_fit_k(tb, cb, c0, bounds, c[-1])[0])
                except Exception:  # noqa: BLE001 - réplica descartada, como en el original
                    pass
            if yb is not None:
                try:
                    popt = _fit_gomp(tb, yb, c0, bounds, y.max())
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
