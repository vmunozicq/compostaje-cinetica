"""Dimensionamiento preliminar de una planta de compostaje (ver docs/modelo_diseno.md).

Une las tres partes del repositorio:

- la cinética de primer orden ajustada a las series de la literatura (k y fracción
  degradable por sustrato),
- el simulador de reactor validado (aire y evaporación por kg de materia orgánica
  degradada, en las condiciones del experimento 2 de Petric y Selimbašić, 2008),
- la estequiometría del mismo artículo (O₂, CO₂ y agua por kg de materia orgánica).

Todas las funciones son puras. Unidades en el nombre de cada campo.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from compostaje import reactor

# Constantes físicas del aire seco a 20 °C y 1 atm.
RHO_AIRE_KG_M3 = 1.204
FRAC_MASICA_O2_AIRE = 0.2315

R2_MIN = 0.5  # mismo filtro de calidad de ajuste que el modelo ML (docs/modelo_ml.md)
PERCENTILES = {"conservador": 0.25, "típico": 0.50, "optimista": 0.75}


@dataclass(frozen=True)
class CineticaSustrato:
    """Resumen de la cinética de primer orden de un sustrato en la literatura."""

    grupo: str
    n_series: int
    k_p25_d: float
    k_p50_d: float
    k_p75_d: float
    frac_degradable: float  # mediana de (C₀ − C∞)/C₀


def cinetica_sustrato(res: pd.DataFrame, grupo: str, r2_min: float = R2_MIN) -> CineticaSustrato:
    """Percentiles de k y fracción degradable mediana de las series de un sustrato.

    `res` es la tabla results/ajustes_por_serie.csv (modo consistente).
    """
    sel = res[(res["Grupo_Residuo"] == grupo) & (res["first_order__r2"] >= r2_min)]
    if sel.empty:
        raise ValueError(f"no hay series con R² ≥ {r2_min} para el sustrato {grupo!r}")
    k = sel["first_order__k"]
    f = (sel["first_order__c0"] - sel["first_order__c_inf"]) / sel["first_order__c0"]
    return CineticaSustrato(grupo, len(sel), float(k.quantile(0.25)), float(k.median()),
                            float(k.quantile(0.75)), float(f.median()))


def tiempo_residencia_d(k_d: float, estabilizacion: float) -> float:
    """Tiempo para degradar la fracción `estabilizacion` de lo degradable (primer orden).

    C(t) = C∞ + (C₀ − C∞)·e^(−k t)  ⇒  t = −ln(1 − X)/k.
    """
    if k_d <= 0:
        raise ValueError("k debe ser positivo")
    if not 0 < estabilizacion < 1:
        raise ValueError("la estabilización debe estar entre 0 y 1 (excluidos)")
    return float(-np.log(1 - estabilizacion) / k_d)


@dataclass(frozen=True)
class RazonesReactor:
    """Razones de operación por kg de materia orgánica degradada, del simulador."""

    aire_m3_por_kg_om: float
    evaporacion_kg_por_kg_om: float
    t_max_c: float


def razones_reactor(op: reactor.Operacion | None = None,
                    par: reactor.Parametros | None = None) -> RazonesReactor:
    """Aire suministrado y agua evaporada por kg de MO degradada en una corrida del simulador.

    Por defecto, el experimento 2 validado (0,18 m³/h, 12,5 kg, 14 días).
    Evaporación = agua inicial − agua final + agua de reacción (y_w · MO degradada).
    """
    op, par = op or reactor.Operacion(), par or reactor.Parametros()
    sim = reactor.simular(op, par)
    om0 = op.masa_kg * (1 - op.humedad) * op.om_base_seca
    om_deg = om0 * sim["conversion_pct"].iloc[-1] / 100
    agua0 = op.masa_kg * op.humedad
    evap = agua0 - sim["m_agua_kg"].iloc[-1] + par.y_w * om_deg
    return RazonesReactor(float(op.caudal_m3_h * op.horas / om_deg), float(evap / om_deg),
                          float(sim["T_C"].max()))


@dataclass(frozen=True)
class EntradaPlanta:
    """Datos de la planta. Densidad y altura no vienen de las fuentes: las da el usuario."""

    residuo_t_d: float  # toneladas húmedas por día
    humedad: float  # fracción, base húmeda
    om_base_seca: float  # fracción de materia orgánica en la materia seca
    estabilizacion: float  # fracción de lo degradable que debe degradarse
    densidad_kg_m3: float  # densidad aparente del material en pila
    altura_pila_m: float

    def __post_init__(self):
        if self.residuo_t_d <= 0:
            raise ValueError("el residuo a tratar debe ser positivo (t/d)")
        for nombre in ("humedad", "om_base_seca", "estabilizacion"):
            v = getattr(self, nombre)
            if not 0 < v < 1 and not (nombre == "om_base_seca" and v == 1):
                raise ValueError(f"{nombre} debe estar entre 0 y 1")
        if self.densidad_kg_m3 <= 0 or self.altura_pila_m <= 0:
            raise ValueError("densidad y altura deben ser positivas")


def dimensionar(ent: EntradaPlanta, cin: CineticaSustrato, raz: RazonesReactor,
                k_usuario_d: float | None = None) -> pd.DataFrame:
    """Dimensiona la planta para los tres escenarios de k (o para un k medido por el usuario).

    Devuelve una fila por escenario. El volumen usa la masa de entrada sin descontar la
    pérdida de masa durante el proceso (lado conservador, supuesto S-D3).
    """
    if k_usuario_d is not None:
        escenarios = {"k medido": k_usuario_d}
    else:
        escenarios = {"conservador (P25)": cin.k_p25_d, "típico (P50)": cin.k_p50_d,
                      "optimista (P75)": cin.k_p75_d}
    entrada_kg_d = ent.residuo_t_d * 1000
    om_entrada_kg_d = entrada_kg_d * (1 - ent.humedad) * ent.om_base_seca
    om_deg_kg_d = om_entrada_kg_d * cin.frac_degradable * ent.estabilizacion
    o2_kg_d = reactor.Parametros().y_o2 * om_deg_kg_d
    aire_esteq_m3_h = o2_kg_d / (FRAC_MASICA_O2_AIRE * RHO_AIRE_KG_M3) / 24
    filas = []
    for nombre, k in escenarios.items():
        t_d = tiempo_residencia_d(k, ent.estabilizacion)
        vol = entrada_kg_d * t_d / ent.densidad_kg_m3
        filas.append({
            "escenario": nombre, "k_d": k, "tiempo_residencia_d": t_d,
            "volumen_m3": vol, "area_pila_m2": vol / ent.altura_pila_m,
            "om_degradada_kg_d": om_deg_kg_d, "o2_kg_d": o2_kg_d,
            "aire_estequiometrico_m3_h": aire_esteq_m3_h,
            "aire_control_termico_m3_h": raz.aire_m3_por_kg_om * om_deg_kg_d / 24,
            "agua_evaporada_kg_d": raz.evaporacion_kg_por_kg_om * om_deg_kg_d,
            "co2_kg_d": reactor.Parametros().y_co2 * om_deg_kg_d,
            "calor_liberado_MJ_d": reactor.Parametros().dh * om_deg_kg_d / 1e6,
        })
    return pd.DataFrame(filas)


def _miles(x: float) -> str:
    return f"{x:,.0f}".replace(",", ".")


def informe_markdown(ent: EntradaPlanta, cin: CineticaSustrato, raz: RazonesReactor,
                     tabla: pd.DataFrame) -> str:
    """Informe descargable con entradas, supuestos y resultados."""
    e = asdict(ent)
    lineas = [
        "# Dimensionamiento preliminar de planta de compostaje",
        "",
        "Generado con github.com/vmunozicq/compostaje-cinetica (docs/modelo_diseno.md).",
        "Es una estimación de prefactibilidad con datos de laboratorio: no reemplaza un",
        "ensayo con el residuo real ni una ingeniería de detalle.",
        "",
        "## Entradas",
        f"- Residuo: {e['residuo_t_d']:g} t/d húmedas, humedad {100 * e['humedad']:.0f} %, "
        f"MO {100 * e['om_base_seca']:.0f} % b.s.",
        f"- Estabilización exigida: {100 * e['estabilizacion']:.0f} % de lo degradable",
        f"- Densidad aparente {e['densidad_kg_m3']:g} kg/m³; altura de pila "
        f"{e['altura_pila_m']:g} m (datos del usuario)",
        f"- Sustrato de referencia: {cin.grupo} ({cin.n_series} series con R² ≥ {R2_MIN}); "
        f"fracción degradable mediana {100 * cin.frac_degradable:.0f} %",
        f"- Aire y evaporación por kg de MO degradada, del simulador validado: "
        f"{raz.aire_m3_por_kg_om:.1f} m³ y {raz.evaporacion_kg_por_kg_om:.2f} kg",
        "",
        "## Resultados",
        "",
        "| Escenario | k (1/d) | Tiempo (d) | Volumen (m³) | Área (m²) | Aire (m³/h) |"
        " Agua evaporada (kg/d) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, r in tabla.iterrows():
        lineas.append(f"| {r.escenario} | {r.k_d:.3f} | {r.tiempo_residencia_d:.0f} | "
                      f"{_miles(r.volumen_m3)} | {_miles(r.area_pila_m2)} | "
                      f"{_miles(r.aire_control_termico_m3_h)} | {_miles(r.agua_evaporada_kg_d)} |")
    r0 = tabla.iloc[0]
    lineas += [
        "",
        f"O₂ consumido {_miles(r0.o2_kg_d)} kg/d (aire estequiométrico "
        f"{_miles(r0.aire_estequiometrico_m3_h)} m³/h); CO₂ {_miles(r0.co2_kg_d)} kg/d; "
        f"calor liberado {_miles(r0.calor_liberado_MJ_d)} MJ/d.",
        "",
        "El aire de diseño es el de control térmico: en el reactor validado se sopla unas "
        f"{raz.aire_m3_por_kg_om / (r0.aire_estequiometrico_m3_h * 24 / r0.om_degradada_kg_d):.0f} "
        "veces el aire estequiométrico, porque el aire retira el calor y el agua.",
    ]
    return "\n".join(lineas) + "\n"
