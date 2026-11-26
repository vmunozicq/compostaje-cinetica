"""Simulador dinámico de un reactor de compostaje (Petric y Selimbašić, 2008).

Especificación, supuestos e interpretaciones en docs/modelo_reactor.md. Las
ecuaciones se citan por su número en el artículo:

    Petric, I., Selimbašić, V. (2008). Development and validation of mathematical
    model for aerobic composting process. Chem. Eng. J. 139(2), 304–317.
    doi:10.1016/j.cej.2007.08.017

Reactor de mezcla completa con tres fases (sólida, líquida y gaseosa): cinética de
orden n con correcciones de temperatura, O₂, humedad y espacio libre de aire;
balances de O₂ y CO₂ disueltos, agua del material, gases (O₂, CO₂, H₂O, N₂) y
energía de las fases gaseosa y sólido-líquida. Unidades internas: kg, kmol, h, K,
Pa, J. Simplificación respecto del artículo: no se modela el NH₃ (requiere el pH
como variable de estado y el propio artículo lo reproduce mal; ver S-R1).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

R = 8314.0  # J kmol⁻¹ K⁻¹
M = {"O2": 32.0, "CO2": 44.0, "H2O": 18.0, "N2": 28.0}  # kg kmol⁻¹
GASES = ("O2", "CO2", "H2O", "N2")
K0 = 273.15


@dataclass(frozen=True)
class Parametros:
    """Constantes de la Tabla 2 y parámetros cinéticos de la Tabla 5 del artículo."""

    # Cinética, Tabla 5 (Ec. 1 y 4): k_T = a·[b^(T−20) − c^(T−60)], T en °C.
    a: float = 0.0000883  # kg^(1−n) h⁻¹
    b: float = 1.0533
    c: float = 1.2247
    n: float = 2.8944
    # Tabla 2. Los c_p de los gases se informan como "J kmol⁻¹ K⁻¹" pero sus valores
    # (29,1; 29,5; …) son J mol⁻¹ K⁻¹: se multiplican por 1000 (interpretación I-1).
    cp_gas: dict = field(default_factory=lambda: {"N2": 29132.0, "O2": 29500.0,
                                                  "CO2": 38154.0, "H2O": 32130.0})
    # c_p de agua, materia orgánica e inorgánica en J kg⁻¹ K⁻¹ (I-1).
    cp_w: float = 4200.0
    cp_om: float = 1320.0
    cp_im: float = 848.0
    kla: float = 1e-4  # kg h⁻¹ Pa⁻¹, igual para O₂, CO₂ y H₂O
    hc: float = 1e5  # J h⁻¹ K⁻¹, transferencia de calor sólido-líquido → gas
    dh: float = 1.54e7  # J kg_OM⁻¹, entalpía de reacción
    ua: float = 1900.0  # J h⁻¹ K⁻¹, pérdidas por las paredes
    k_o2: float = 0.07  # kg O₂ m⁻³, constante de saturación de O₂
    y_o2: float = 0.705  # kg O₂ / kg OM
    y_co2: float = 1.276  # kg CO₂ / kg OM
    y_w: float = 0.360  # kg H₂O / kg OM
    # Espacio libre de aire (Ec. 8–10). C no tiene valor en el artículo (rango 0,15–0,4):
    # supuesto S-R3.
    c_bulk: float = 0.25
    g_v: float = 1.0
    g_f: float = 2.5
    # Calor latente del agua para la entalpía del vapor (supuesto S-R2).
    lambda_w: float = 2.44e6  # J kg⁻¹
    # Qué factores de la Ec. 2 se aplican (hallazgo H-1 de docs/modelo_reactor.md):
    # "articulo" usa solo k_T, que es lo que reproduce las curvas publicadas con los
    # parámetros de la Tabla 5; "todas" aplica también k_O2, k_H2O y k_FAS.
    correcciones: str = "articulo"


@dataclass(frozen=True)
class Operacion:
    """Condiciones del experimento 2 del artículo (Figs. 3–6), por defecto."""

    masa_kg: float = 12.5
    humedad: float = 0.6911  # fracción, base húmeda
    om_base_seca: float = 0.8022  # fracción
    caudal_m3_h: float = 0.18
    t0_c: float = 28.0  # supuesto S-R4 (leído de la Fig. 3)
    ta_c: float = 28.0  # supuesto S-R4
    v_reactor_m3: float = 0.032
    p_total: float = 101325.0
    horas: float = 336.0

    def __post_init__(self):
        if not 0 < self.humedad < 1 or not 0 < self.om_base_seca <= 1:
            raise ValueError("humedad y fracción orgánica deben estar entre 0 y 1")
        if self.masa_kg <= 0 or self.caudal_m3_h < 0:
            raise ValueError("masa positiva y caudal no negativo")


def p_sat(t_k):
    """Presión de vapor del agua (Pa), Ec. 26."""
    return 10 ** (22.443 - 2795.0 / t_k - 1.6798 * np.log(t_k))


def henry(gas: str, t_k):
    """Constante de Henry (Pa), Ec. 15–16. La Ec. 16 da atm para el CO₂ (I-2)."""
    if gas == "O2":
        return 101325.0 * np.exp(66.7354 - 8747.55 / t_k - 24.4526 * np.log(t_k / 100.0))
    if gas == "CO2":
        return 101325.0 * (-11418.84 + 43.8658 * t_k)
    raise ValueError(gas)


def _fas(s_m, v_s, par: Parametros):
    """Espacio libre de aire (Ec. 8–10), con δ_w = 1 en las mismas unidades que C."""
    delta_m = par.c_bulk / s_m
    g_s = 1.0 / (v_s / par.g_v + (1 - v_s) / par.g_f)
    return 1 - delta_m * s_m / g_s - delta_m * (1 - s_m)


def correcciones(t_c, o2_kg_m3, s_m, v_s, par: Parametros) -> dict:
    """Factores de la Ec. 2: k_T (Ec. 4), k_O₂ (Ec. 5), k_H₂O (Ec. 6), k_FAS (Ec. 7)."""
    k_t = max(par.a * (par.b ** (t_c - 20) - par.c ** (t_c - 60)), 0.0)
    k_o2 = max(o2_kg_m3, 0.0) / (par.k_o2 + max(o2_kg_m3, 0.0))
    k_h2o = 1.0 / (np.exp(-17.684 * (1 - s_m) + 7.0622) + 1)
    k_fas = 1.0 / (np.exp(-23.675 * _fas(s_m, v_s, par) + 3.4945) + 1)
    return {"k_T": k_t, "k_O2": k_o2, "k_H2O": k_h2o, "k_FAS": k_fas}


def _setup(op: Operacion, par: Parametros):
    seca = op.masa_kg * (1 - op.humedad)
    m_om, m_im, m_w = seca * op.om_base_seca, seca * (1 - op.om_base_seca), op.masa_kg * op.humedad
    s_m = 1 - op.humedad
    rho_kg_m3 = 1000.0 * par.c_bulk / s_m  # δ_m (Ec. 10) en t m⁻³ → kg m⁻³
    v_g = op.v_reactor_m3 - op.masa_kg / rho_kg_m3
    if v_g <= 0:
        raise ValueError("el material no cabe en el reactor con esa densidad")
    psi0 = op.ta_c + K0  # el aire entra a temperatura ambiente
    ps0 = p_sat(psi0)
    q = op.caudal_m3_h
    f0 = {"O2": 0.21 * (op.p_total - ps0) * q / (R * psi0),  # Ec. 28
          "N2": 0.79 * (op.p_total - ps0) * q / (R * psi0),  # Ec. 29
          "H2O": ps0 * q / (R * psi0),  # Ec. 30
          "CO2": 0.0}
    t0 = op.t0_c + K0
    n_tot = op.p_total * v_g / (R * t0)
    ps_t0 = p_sat(t0)
    n0 = {"O2": 0.21 * (n_tot - ps_t0 * v_g / (R * t0)),
          "N2": 0.79 * (n_tot - ps_t0 * v_g / (R * t0)),
          "H2O": ps_t0 * v_g / (R * t0), "CO2": 1e-9}
    # Gases disueltos en equilibrio con el aire inicial (Ec. 14 con R^T = 0).
    kmol_w = m_w / M["H2O"]
    x_o2 = n0["O2"] * R * t0 / v_g / henry("O2", t0)
    m_o2d = x_o2 * kmol_w * M["O2"]
    y0 = np.array([m_om, m_o2d, 0.0, m_w, n0["O2"], n0["CO2"], n0["H2O"], n0["N2"], t0, t0])
    return y0, {"m_im": m_im, "v_g": v_g, "psi0": psi0, "f0": f0, "q": q}


def _rhs(t, y, op: Operacion, par: Parametros, ctx: dict, diag: bool = False):
    m_om, m_o2d, m_co2d, m_w, n_o2, n_co2, n_h2o, n_n2, psi, temp = y
    m_om, m_w = max(m_om, 1e-9), max(m_w, 1e-9)
    n = {"O2": max(n_o2, 0.0), "CO2": max(n_co2, 0.0), "H2O": max(n_h2o, 0.0),
         "N2": max(n_n2, 0.0)}
    v_g, m_im = ctx["v_g"], ctx["m_im"]
    p = {g: n[g] * R * psi / v_g for g in GASES}  # Ec. 22
    s_m = (m_om + m_im) / (m_om + m_im + m_w)
    v_s = m_om / (m_om + m_im)
    corr = correcciones(temp - K0, n["O2"] * M["O2"] / v_g, s_m, v_s, par)
    if par.correcciones == "todas":
        k = corr["k_T"] * corr["k_O2"] * corr["k_H2O"] * corr["k_FAS"]  # Ec. 2
    elif par.correcciones == "articulo":
        k = corr["k_T"]
    else:
        raise ValueError("correcciones debe ser 'articulo' o 'todas'")
    r = k * m_om ** par.n  # velocidad de degradación, kg h⁻¹ (Ec. 1)
    dm_om = -r

    # Transferencia líquido → gas, kg h⁻¹ (Ec. 14 y 25); positiva hacia el gas.
    kmol_w = m_w / M["H2O"]
    x_o2 = max(m_o2d, 0.0) / M["O2"] / kmol_w
    x_co2 = max(m_co2d, 0.0) / M["CO2"] / kmol_w
    rt = {"O2": par.kla * (henry("O2", temp) * x_o2 - p["O2"]),
          "CO2": par.kla * (henry("CO2", temp) * x_co2 - p["CO2"]),
          "H2O": par.kla * (p_sat(psi) - p["H2O"]),
          "N2": 0.0}
    # Generación por reacción (Ec. 13 y 24), con el signo físico (I-3): se consume O₂,
    # se producen CO₂ y agua.
    dm_o2d = -par.y_o2 * r - rt["O2"]
    dm_co2d = par.y_co2 * r - rt["CO2"]
    dm_w = par.y_w * r - rt["H2O"]
    rbar = {g: rt[g] / M[g] for g in GASES}  # kmol h⁻¹

    # Temperatura del gas (Ec. 32–33).
    cp = par.cp_gas
    f0 = ctx["f0"]
    qc = par.hc * (temp - psi)
    dpsi = ((qc + (ctx["psi0"] - psi) * sum(cp[g] * f0[g] for g in GASES)
             + (temp - psi) * sum(cp[g] * max(0.0, rbar[g]) for g in GASES))
            / max(sum(n[g] * cp[g] for g in GASES), 1e-12))

    # Gases en la fase gaseosa (Ec. 27 y 31).
    f_t0 = sum(f0.values())
    flujo = R * psi / (op.p_total * v_g) * (f_t0 + sum(rbar.values()))
    dn = {g: f0[g] + rbar[g] - (n[g] * flujo + n[g] / psi * dpsi) for g in GASES}

    # Temperatura sólido-líquida (Ec. 34, 36–38). Entalpía molar del gas: c_p·(T − 273,15),
    # más el calor latente para el vapor de agua (supuesto S-R2).
    def h(g, tk):
        lat = par.lambda_w * M["H2O"] if g == "H2O" else 0.0
        return cp[g] * (tk - K0) + lat

    q_cw = par.ua * (op.ta_c + K0 - temp)
    q_g = par.dh * r
    gas_term = sum(max(0.0, rbar[g]) * h(g, temp) + min(0.0, rbar[g]) * h(g, psi)
                   for g in GASES)
    dtemp = ((q_cw + q_g - qc - gas_term - (temp - K0) * (par.cp_w * dm_w + par.cp_om * dm_om))
             / (par.cp_w * m_w + par.cp_om * m_om + par.cp_im * m_im))
    dy = [dm_om, dm_o2d, dm_co2d, dm_w, dn["O2"], dn["CO2"], dn["H2O"], dn["N2"], dpsi, dtemp]
    if diag:
        return dy, {"r": r, "k": k, **corr, "Q_G": q_g, "Q_cw": q_cw}
    return dy


def simular(op: Operacion | None = None, par: Parametros | None = None,
            paso_h: float = 1.0) -> pd.DataFrame:
    """Integra el modelo y devuelve una fila por hora (por defecto)."""
    op, par = op or Operacion(), par or Parametros()
    y0, ctx = _setup(op, par)
    t_eval = np.arange(0.0, op.horas + 1e-9, paso_h)
    sol = solve_ivp(_rhs, (0.0, op.horas), y0, method="LSODA", t_eval=t_eval,
                    args=(op, par, ctx), rtol=1e-6, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"la integración falló: {sol.message}")
    m_om, _, _, m_w, n_o2, n_co2, _, n_n2, psi, temp = sol.y
    seco = n_o2 + n_co2 + n_n2
    out = pd.DataFrame({
        "t_h": sol.t, "T_C": temp - K0, "T_gas_C": psi - K0, "m_OM_kg": m_om,
        "conversion_pct": 100 * (1 - m_om / y0[0]), "m_agua_kg": m_w,
        "humedad_pct": 100 * m_w / (m_w + m_om + ctx["m_im"]),
        # Gas de salida en base seca (medido tras el condensador, I-4).
        "O2_pct": 100 * n_o2 / seco, "CO2_pct": 100 * n_co2 / seco,
    })
    return out


def variar(campo: str, valores, op: Operacion | None = None,
           par: Parametros | None = None) -> pd.DataFrame:
    """Simula variando un campo de Operacion (p. ej. caudal o humedad), como Figs. 8–11."""
    op = op or Operacion()
    filas = []
    for v in valores:
        s = simular(replace(op, **{campo: v}), par)
        filas.append({campo: v, "T_max_C": s["T_C"].max(),
                      "t_T_max_h": s.loc[s["T_C"].idxmax(), "t_h"],
                      "conversion_final_pct": s["conversion_pct"].iloc[-1],
                      "O2_min_pct": s["O2_pct"].min(),
                      "horas_sobre_55C": float((s["T_C"] >= 55).sum())})
    return pd.DataFrame(filas)


def cargar_experimento(path) -> pd.DataFrame:
    """Datos digitalizados del experimento 2 (data/reference/petric2008_exp2_digitalizado.csv)."""
    return pd.read_csv(path, comment="#")


def validar(sim: pd.DataFrame, exp: pd.DataFrame) -> pd.DataFrame:
    """Error del simulador frente a cada variable medida: RMSE, sesgo medio y máximo."""
    filas = []
    for var, g in exp.groupby("variable"):
        pred = np.interp(g["t_h"], sim["t_h"], sim[var])
        err = pred - g["valor"].to_numpy()
        filas.append({"variable": var, "n": len(g), "RMSE": float(np.sqrt(np.mean(err ** 2))),
                      "sesgo": float(err.mean()), "error_max": float(np.abs(err).max())})
    return pd.DataFrame(filas)
