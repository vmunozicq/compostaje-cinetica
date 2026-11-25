"""Métricas de bondad de ajuste y criterios de información (Anexo B de la memoria)."""

from __future__ import annotations

import numpy as np


def sse(y, yhat) -> float:
    r = np.asarray(y, float) - np.asarray(yhat, float)
    return float(np.sum(r**2))


def sst(y) -> float:
    y = np.asarray(y, float)
    return float(np.sum((y - y.mean()) ** 2))


def r2(y, yhat) -> float:
    """R² = 1 − SSE/SST. Devuelve NaN si la serie es plana (SST < 1e-5)."""
    s = sst(y)
    if s < 1e-5:
        return float("nan")
    return 1.0 - sse(y, yhat) / s


def r2_adj(r2_value: float, n: int, p: int) -> float:
    if n - p - 1 <= 0:
        return float("nan")
    return 1.0 - (1.0 - r2_value) * (n - 1) / (n - p - 1)


def rmse(y, yhat) -> float:
    return float(np.sqrt(np.mean((np.asarray(y, float) - np.asarray(yhat, float)) ** 2)))


def nrmse_pct(y, yhat) -> float:
    """RMSE normalizado por la media observada, en % (definición del Anexo B)."""
    return 100.0 * rmse(y, yhat) / float(np.mean(y))


def information_criteria(sse_value: float, n: int, p: int) -> tuple[float, float, float]:
    """AIC = n·ln(SSE/n) + 2p; AICc = AIC + 2p(p+1)/(n−p−1); BIC = n·ln(SSE/n) + p·ln(n).

    Devuelve (inf, inf, inf) si n ≤ p + 1 (AICc indefinido), igual que el script original.
    """
    if not np.isfinite(sse_value) or n <= p + 1:
        return (np.inf, np.inf, np.inf)
    s = max(sse_value, 1e-10)
    base = n * np.log(s / n)
    aic = base + 2 * p
    aicc = aic + 2 * p * (p + 1) / (n - p - 1)
    bic = base + p * np.log(n)
    return (float(aic), float(aicc), float(bic))
