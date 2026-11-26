"""Modelos cinéticos de la memoria (formas usadas en el ajuste final).

Notación: C(t) es la concentración de TOC u OM (%), C0 = primer dato medido
(fijo, no se ajusta) y C_inf es la fracción no degradable (asíntota).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import odeint

R_GAS = 8.314  # J mol^-1 K^-1


def first_order(t, k, c_inf, c0):
    """C(t) = C_inf + (C0 − C_inf)·exp(−k·t)."""
    t = np.asarray(t, float)
    return c_inf + (c0 - c_inf) * np.exp(-k * t)


def biphasic(t, k1, k2, f, c_inf, c0):
    """C(t) = C_inf + (C0 − C_inf)·[f·exp(−k1·t) + (1 − f)·exp(−k2·t)]."""
    t = np.asarray(t, float)
    return c_inf + (c0 - c_inf) * (f * np.exp(-k1 * t) + (1 - f) * np.exp(-k2 * t))


def gompertz_degraded(t, p, rm, lam):
    """Degradación acumulada Y(t) = C0 − C(t) (Gompertz modificado, Zwietering 1990).

    Y(t) = P·exp(−exp(Rm·e/P·(λ − t) + 1))
    """
    t = np.asarray(t, float)
    p = max(p, 1e-12)
    arg = np.clip(rm * np.e / p * (lam - t) + 1.0, -700, 700)
    return p * np.exp(-np.exp(arg))


def gompertz(t, p, rm, lam, c0):
    """C(t) = C0 − Y(t)."""
    return c0 - gompertz_degraded(t, p, rm, lam)


def k_arrhenius(temp_c, a, ea):
    """k(T) = A·exp(−Ea / (R·T)), T en K."""
    return a * np.exp(-ea / (R_GAS * (np.asarray(temp_c, float) + 273.15)))


def arrhenius(t_eval, a, ea, c_inf, c0, t_temp, temp_c):
    """Resuelve dC/dt = −k(T(t))·(C − C_inf), C(0) = C0.

    T(t) se interpola linealmente entre las mediciones de temperatura
    (``np.interp``: fuera del rango se mantiene el valor extremo).
    """
    t_temp = np.asarray(t_temp, float)
    temp_c = np.asarray(temp_c, float)

    def rhs(c, t):
        if c[0] <= c_inf:
            return [0.0]
        return [-k_arrhenius(np.interp(t, t_temp, temp_c), a, ea) * (c[0] - c_inf)]

    t_eval = np.asarray(t_eval, float)
    # odeint integra desde el primer tiempo pedido; se antepone t0 del dato.
    t0 = min(t_temp[0], t_eval[0])
    grid = np.concatenate([[t0], t_eval]) if t_eval[0] > t0 else t_eval
    sol = odeint(rhs, [c0], grid).ravel()
    return sol[-len(t_eval):]


def k_arrhenius_tref(temp_c, k_ref, ea, t_ref_c):
    """k(T) = k_ref·exp[−Ea/R·(1/T − 1/T_ref)], T y T_ref en °C (se convierten a K).

    Reparametrización estándar de Arrhenius (docs/supuestos.md, L1): k_ref es la
    constante a T_ref y queda mucho menos correlacionada con Ea que A.
    """
    t_k = np.asarray(temp_c, float) + 273.15
    return k_ref * np.exp(-ea / R_GAS * (1.0 / t_k - 1.0 / (t_ref_c + 273.15)))


def a_from_kref(k_ref, ea, t_ref_c):
    """Factor preexponencial equivalente: A = k_ref·exp(Ea / (R·T_ref))."""
    return k_ref * np.exp(ea / (R_GAS * (t_ref_c + 273.15)))


def arrhenius_tref(t_eval, k_ref, ea, c_inf, c0, t_temp, temp_c, t_ref_c):
    """Mismo ODE que :func:`arrhenius`, con k(T) en la forma de temperatura de referencia."""
    return arrhenius(t_eval, a_from_kref(k_ref, ea, t_ref_c), ea, c_inf, c0, t_temp, temp_c)


# Número de parámetros ajustables por modelo (para AIC/AICc/BIC), como en la memoria.
N_PARAMS = {"first_order": 2, "biphasic": 4, "gompertz": 3, "arrhenius": 3}
MODEL_LABELS = {
    "first_order": "Primer orden",
    "biphasic": "Bifásico",
    "gompertz": "Gompertz modificado",
    "arrhenius": "Arrhenius",
}
