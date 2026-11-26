"""Ajuste por regresión no lineal de los cuatro modelos a una serie.

Valores iniciales y cotas tomados de los scripts que generaron los resultados de
la memoria (``modelacion_masiva.py``). Ver ``docs/supuestos.md``.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit, least_squares, minimize

from . import metrics
from .models import (
    N_PARAMS,
    a_from_kref,
    arrhenius,
    arrhenius_tref,
    biphasic,
    first_order,
    gompertz,
    gompertz_degraded,
)

# Cotas del ajuste principal (modelacion_masiva.py).
BOUNDS = {
    "k_max": 10.0,          # primer orden y bifásico, d^-1
    "rm_max": 100.0,        # Gompertz, %·d^-1
    "lam_max": 50.0,        # Gompertz, d
    "a_range": (1.0, 1e9),  # Arrhenius, d^-1
    "ea_range": (1000.0, 150000.0),  # Arrhenius, J/mol
}


@dataclass
class FitResult:
    model: str
    params: dict
    success: bool
    n: int
    p: int
    sse: float = np.nan
    r2: float = np.nan
    r2_adj: float = np.nan
    rmse: float = np.nan
    nrmse_pct: float = np.nan
    aic: float = np.inf
    aicc: float = np.inf
    bic: float = np.inf
    message: str = ""
    predict: object = field(default=None, repr=False)

    def as_row(self) -> dict:
        row = {f"{k}": v for k, v in self.params.items()}
        row.update(dict(r2=self.r2, r2_adj=self.r2_adj, rmse=self.rmse,
                        nrmse_pct=self.nrmse_pct, aic=self.aic, aicc=self.aicc, bic=self.bic))
        return {f"{self.model}__{k}": v for k, v in row.items()}


def _finish(model, params, t, c, predict) -> FitResult:
    p = N_PARAMS[model]
    yhat = predict(t)
    s = metrics.sse(c, yhat)
    r2 = metrics.r2(c, yhat)
    res = FitResult(model=model, params=params, success=True, n=len(t), p=p, sse=s, r2=r2,
                    r2_adj=metrics.r2_adj(r2, len(t), p), rmse=metrics.rmse(c, yhat),
                    nrmse_pct=metrics.nrmse_pct(c, yhat), predict=predict)
    # Igual que el script de criterios: un ajuste con R² ≤ 0 no compite.
    if np.isfinite(r2) and r2 > 0:
        res.aic, res.aicc, res.bic = metrics.information_criteria(s, len(t), p)
    return res


def _failed(model, n, msg) -> FitResult:
    return FitResult(model=model, params={}, success=False, n=n, p=N_PARAMS[model], message=msg)


def fit_first_order(t, c, bounds=BOUNDS) -> FitResult:
    c0 = c[0]
    try:
        popt, _ = curve_fit(lambda tt, k, ci: first_order(tt, k, ci, c0), t, c,
                            p0=[0.1, c[-1]], bounds=([0, 0], [bounds["k_max"], c0]))
    except Exception as e:  # noqa: BLE001 - se registra el motivo
        return _failed("first_order", len(t), str(e))
    k, ci = popt
    return _finish("first_order", {"k": k, "c_inf": ci, "c0": c0}, t, c,
                   lambda tt: first_order(tt, k, ci, c0))


def fit_biphasic(t, c, bounds=BOUNDS) -> FitResult:
    c0 = c[0]
    try:
        popt, _ = curve_fit(lambda tt, k1, k2, f, ci: biphasic(tt, k1, k2, f, ci, c0), t, c,
                            p0=[0.5, 0.05, 0.5, c[-1]],
                            bounds=([0, 0, 0, 0], [bounds["k_max"], bounds["k_max"], 1, c0]))
    except Exception as e:  # noqa: BLE001
        return _failed("biphasic", len(t), str(e))
    k1, k2, f, ci = popt
    pred = lambda tt: biphasic(tt, k1, k2, f, ci, c0)  # noqa: E731
    # Etiquetado: la fracción "rápida" es la de mayor k (el modelo es simétrico).
    if k2 > k1:
        k1, k2, f = k2, k1, 1 - f
    return _finish("biphasic", {"k1": k1, "k2": k2, "f": f, "c_inf": ci, "c0": c0}, t, c, pred)


def fit_gompertz(t, c, bounds=BOUNDS) -> FitResult:
    c0 = c[0]
    y = c0 - c  # degradación acumulada
    try:
        popt, _ = curve_fit(gompertz_degraded, t, y, p0=[y.max(), 1, 1],
                            bounds=([0, 0, 0], [c0, bounds["rm_max"], bounds["lam_max"]]),
                            maxfev=5000)
    except Exception as e:  # noqa: BLE001
        return _failed("gompertz", len(t), str(e))
    p, rm, lam = popt
    return _finish("gompertz", {"P": p, "Rm": rm, "lambda": lam, "c0": c0}, t, c,
                   lambda tt: gompertz(tt, p, rm, lam, c0))


def fit_arrhenius(t, c, temp, bounds=BOUNDS) -> FitResult:
    c0 = c[0]
    if temp is None or np.isnan(temp).any():
        return _failed("arrhenius", len(t), "sin perfil de temperatura")

    def obj(x):
        return np.sum((c - arrhenius(t, x[0], x[1], x[2], c0, t, temp)) ** 2)

    try:
        res = minimize(obj, [1e4, 30000, c[-1]], method="L-BFGS-B",
                       bounds=[bounds["a_range"], bounds["ea_range"], (0, c0)])
    except Exception as e:  # noqa: BLE001
        return _failed("arrhenius", len(t), str(e))
    a, ea, ci = res.x
    t_temp, temp_c = np.array(t), np.array(temp)
    return _finish("arrhenius", {"A": a, "Ea": ea, "c_inf": ci, "c0": c0}, t, c,
                   lambda tt: arrhenius(tt, a, ea, ci, c0, t_temp, temp_c))


# Valores iniciales de Ea (J/mol) para el ajuste multiarranque de Arrhenius con T_ref.
EA_STARTS = (20000.0, 60000.0, 100000.0)


def fit_arrhenius_tref(t, c, temp, bounds=BOUNDS) -> FitResult:
    """Arrhenius con temperatura de referencia: k(T) = k_ref·exp[−Ea/R·(1/T − 1/T_ref)].

    T_ref es la temperatura media medida en la serie (°C), de modo que k_ref es la
    constante a la temperatura típica del ensayo. Mismo ODE, mismos p = 3 parámetros
    (k_ref, Ea, C_inf) y mismo rango de Ea que la forma original. Se ajusta con
    mínimos cuadrados (TRF con cotas) desde varios Ea iniciales y se guarda el de
    menor SSE. Se informa también el A equivalente.
    """
    c0 = c[0]
    if temp is None or np.isnan(temp).any():
        return _failed("arrhenius", len(t), "sin perfil de temperatura")
    t_temp, temp_c = np.array(t), np.array(temp)
    t_ref = float(np.mean(temp_c))

    def resid(x):
        return c - arrhenius_tref(t, x[0], x[1], x[2], c0, t_temp, temp_c, t_ref)

    lo = [0.0, bounds["ea_range"][0], 0.0]
    hi = [bounds["k_max"], bounds["ea_range"][1], c0]
    best = None
    for ea0 in EA_STARTS:
        try:
            r = least_squares(resid, [0.1, ea0, min(c[-1], c0)], bounds=(lo, hi),
                              x_scale=[0.1, 1e4, max(c0, 1.0)], max_nfev=300)
        except Exception:  # noqa: BLE001 - se prueba el siguiente arranque
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        return _failed("arrhenius", len(t), "no convergió desde ningún arranque")
    k_ref, ea, ci = best.x
    return _finish("arrhenius", {"k_ref": k_ref, "Ea": ea, "c_inf": ci, "T_ref": t_ref,
                                 "A": a_from_kref(k_ref, ea, t_ref), "c0": c0}, t, c,
                   lambda tt: arrhenius_tref(tt, k_ref, ea, ci, c0, t_temp, temp_c, t_ref))


FITTERS = {
    "first_order": fit_first_order,
    "biphasic": fit_biphasic,
    "gompertz": fit_gompertz,
}


MODOS = ("consistente", "memoria")


def _shifted(predict, t0):
    return lambda tt: predict(np.asarray(tt, float) - t0)


def fit_all(t, c, temp=None, bounds=BOUNDS, modo: str = "consistente") -> dict[str, FitResult]:
    """Ajusta los cuatro modelos a una serie.

    modo="memoria" reproduce los scripts de la memoria: el modelo impone C(0) = C0 aunque
    la serie parta en t > 0, y Arrhenius usa la forma A·exp(−Ea/RT) con L-BFGS-B.
    modo="consistente" (por defecto) mide el tiempo desde el primer dato (t − t0, L2) y
    ajusta Arrhenius con temperatura de referencia (L1). Ver docs/supuestos.md.
    Las funciones ``predict`` reciben siempre el tiempo original de la serie.
    """
    if modo not in MODOS:
        raise ValueError(f"modo debe ser uno de {MODOS}, no {modo!r}")
    t = np.asarray(t, float)
    c = np.asarray(c, float)
    if t.shape != c.shape:
        raise ValueError("t y C deben tener el mismo largo")
    t0 = float(t[0]) if modo == "consistente" else 0.0
    ts = t - t0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (OptimizeWarning, RuntimeWarning))
        out = {name: f(ts, c, bounds) for name, f in FITTERS.items()}
        temp_arr = None if temp is None else np.asarray(temp, float)
        fit_arr = fit_arrhenius_tref if modo == "consistente" else fit_arrhenius
        out["arrhenius"] = fit_arr(ts, c, temp_arr, bounds)
    if t0:
        for r in out.values():
            if r.predict is not None:
                r.predict = _shifted(r.predict, t0)
    return out


def best_model(fits: dict[str, FitResult], criterion: str = "aicc") -> str | None:
    vals = {m: getattr(r, criterion) for m, r in fits.items()}
    finite = {m: v for m, v in vals.items() if np.isfinite(v)}
    if not finite:
        return None
    return min(finite, key=finite.get)
