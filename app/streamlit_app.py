"""Interfaz Streamlit: explorar el ajuste de los 4 modelos cinéticos de compostaje."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from compostaje import bootstrap, data, fitting, ml, pipeline, reactor  # noqa: E402
from compostaje.models import MODEL_LABELS  # noqa: E402
from compostaje.plots import GROUP_COLORS, GROUP_LABELS, MODEL_COLORS  # noqa: E402

DASH = {"first_order": "dash", "biphasic": "dashdot", "gompertz": "solid", "arrhenius": "dot"}
SHORT_ES = {"1st": "Primer orden", "Gomp": "Gompertz", "Bi": "Bifásico", "Arr": "Arrhenius",
            "Ninguno": "Ninguno"}

st.set_page_config(page_title="Cinética de compostaje", page_icon="🌱", layout="wide")


@st.cache_data(show_spinner="Cargando dataset…")
def load_valid() -> pd.DataFrame:
    return data.build_validated_dataset(data.load_raw(data.DEFAULT_RAW_PATH))


@st.cache_data
def load_results(modo: str = "consistente") -> pd.DataFrame:
    nombre = "ajustes_por_serie.csv" if modo == "consistente" else "ajustes_modo_memoria.csv"
    return pd.read_csv(ROOT / "results" / nombre, index_col="ID")


@st.cache_resource(show_spinner="Ajustando modelos…", max_entries=64)
def fit_series(t: tuple, c: tuple, temp: tuple | None, modo: str = "consistente"):
    # cache_resource: los resultados incluyen funciones de predicción (no serializables).
    return fitting.fit_all(np.array(t), np.array(c), None if temp is None else np.array(temp),
                           modo=modo)


@st.cache_data(show_spinner="Bootstrap en curso…")
def run_bootstrap(t: tuple, c: tuple, n_boot: int, metodo: str, modo: str):
    out = bootstrap.bootstrap_series(np.array(t), np.array(c), n_boot=n_boot, seed=42,
                                     metodo=metodo, modo=modo)
    return {k: {q: v[q] for q in ("mean", "lo", "hi", "n_ok")} for k, v in out.items()}


@st.cache_resource(show_spinner="Entrenando el modelo de ML…")
def ml_model():
    tab = pd.read_csv(ROOT / "results" / "ml_tabla.csv", index_col="ID")
    model, dom = ml.fit_final(tab)
    return tab, model, dom


@st.cache_data(show_spinner="Integrando el balance de masa y energía…")
def simular_reactor(caudal: float, humedad: float, masa: float, t0: float, ta: float,
                    correcciones: str) -> pd.DataFrame:
    op = reactor.Operacion(masa_kg=masa, humedad=humedad, caudal_m3_h=caudal, t0_c=t0, ta_c=ta)
    return reactor.simular(op, reactor.Parametros(correcciones=correcciones))


def fit_figure(t, c, fits, ylabel, title):
    fig = go.Figure()
    t_s = np.linspace(min(t), max(t), 200)
    for m in ("first_order", "biphasic", "gompertz", "arrhenius"):
        r = fits[m]
        if not r.success or not np.isfinite(r.r2):
            continue
        fig.add_trace(go.Scatter(
            x=t_s, y=r.predict(t_s), mode="lines", name=f"{MODEL_LABELS[m]} (R² {r.r2:.2f})",
            line=dict(color=MODEL_COLORS[m], width=2, dash=DASH[m]),
            hovertemplate="t = %{x:.1f} d<br>C = %{y:.2f} %<extra>" + MODEL_LABELS[m] + "</extra>"))
    fig.add_trace(go.Scatter(
        x=t, y=c, mode="markers", name="Datos experimentales",
        marker=dict(color="#191919", size=9, line=dict(color="white", width=1.5)),
        hovertemplate="t = %{x:.1f} d<br>C = %{y:.2f} %<extra>Dato</extra>"))
    fig.update_layout(title=title, xaxis_title="Tiempo (d)", yaxis_title=ylabel,
                      template="plotly_white", height=430, hovermode="closest",
                      legend=dict(orientation="h", y=-0.2), margin=dict(t=50, b=10))
    return fig


def metrics_table(fits) -> pd.DataFrame:
    best = fitting.best_model(fits)
    rows = []
    for m, r in fits.items():
        if not r.success:
            rows.append({"Modelo": MODEL_LABELS[m], "Parámetros": "no convergió"})
            continue
        params = ", ".join(f"{k} = {v:.4g}" for k, v in r.params.items() if k != "c0")
        rows.append({"Modelo": MODEL_LABELS[m] + ("  ★" if m == best else ""),
                     "Parámetros": params, "R²": r.r2, "RMSE (%)": r.rmse,
                     "AICc": r.aicc if np.isfinite(r.aicc) else None})
    return pd.DataFrame(rows)


st.title("Modelación cinética del compostaje aeróbico")
st.caption("Primer orden · Bifásico · Gompertz modificado · Arrhenius, con selección por AICc "
           "e incertidumbre bootstrap. Basado en mi memoria de título (Ing. Civil Químico, UdeC).")

tab1, tab2, tab3, tab_ml, tab_r, tab4 = st.tabs([
    "Explorar una serie", "Resumen de 159 series", "Ajustar mis datos", "Predecir k (ML)",
    "Simular reactor", "Acerca de"])

with tab1:
    valid = load_valid()
    res = load_results()
    col_a, col_b = st.columns([1, 3])
    with col_a:
        grupos = [g for g in GROUP_LABELS if g in set(res["Grupo_Residuo"])]
        grupo = st.selectbox("Sustrato", grupos, format_func=GROUP_LABELS.get)
        ids = res.index[res["Grupo_Residuo"] == grupo].tolist()
        default = {"Kitchen Waste": "25_1", "Sludge": "461_1", "Manure": "152_1"}.get(grupo)
        uid = st.selectbox("Serie (ID)", ids, index=ids.index(default) if default in ids else 0)
        st.caption("Casos de estudio de la memoria: 25_1 (cocina), 461_1 (lodos), 152_1 "
                   "(estiércol).")
        modo = st.radio("Cálculo", ["consistente", "memoria"], horizontal=True,
                        format_func={"consistente": "Versión corregida",
                                     "memoria": "Como en la memoria"}.get,
                        help="Corregida: tiempo desde el primer dato (t − t₀) y Arrhenius con "
                             "temperatura de referencia. Memoria: reproduce los scripts "
                             "originales. Ver docs/supuestos.md.")
        metodo = st.radio("Bootstrap", ["residuos", "pares"], horizontal=True,
                          format_func={"residuos": "De residuos", "pares": "De pares"}.get,
                          help="Residuos: más estable con pocos puntos. Pares: método de la "
                               "memoria.")
        n_boot = st.select_slider("Réplicas bootstrap", [100, 250, 500, 1000], value=250)
        do_boot = st.button("Calcular IC 95 % (bootstrap)")
    s = data.extract_series(valid[valid["Unique_ID"] == uid])
    fits = fit_series(tuple(s.t), tuple(s.c), tuple(s.temp), modo)
    ylabel = "TOC (%)" if s.variable == "TOC_percent" else "Materia orgánica (%)"
    with col_b:
        st.plotly_chart(fit_figure(s.t, s.c, fits, ylabel,
                                   f"Serie {uid} · {s.waste_type} · n = {s.n}"),
                        use_container_width=True)
    st.dataframe(metrics_table(fits), hide_index=True, use_container_width=True,
                 column_config={"R²": st.column_config.NumberColumn(format="%.3f"),
                                "RMSE (%)": st.column_config.NumberColumn(format="%.3f"),
                                "AICc": st.column_config.NumberColumn(format="%.2f")})
    st.caption("★ = menor AICc (mejor equilibrio entre ajuste y n.º de parámetros). "
               "Con n ≤ p + 1 el AICc no está definido y el modelo no compite.")
    if do_boot:
        b = run_bootstrap(tuple(s.t), tuple(s.c), n_boot, metodo, modo)
        c1, c2, c3 = st.columns(3)
        for col, key, label, unit in ((c1, "k", "k primer orden", "d⁻¹"),
                                      (c2, "Rm", "Rm Gompertz", "%/d"),
                                      (c3, "lambda", "λ Gompertz", "d")):
            v = b[key]
            col.metric(label, f"{v['mean']:.3g} {unit}",
                       help=f"IC 95 %: {v['lo']:.3g} – {v['hi']:.3g} "
                            f"({v['n_ok']} réplicas convergidas)")
            col.caption(f"IC 95 %: {v['lo']:.3g} – {v['hi']:.3g}")
        st.info("Para diseño conservador, la memoria recomienda usar el límite inferior del "
                "IC de k, no la media.")

with tab2:
    res = load_results()
    counts = pipeline.winner_counts(res, load_results("memoria"))
    st.subheader("¿Qué modelo gana según AICc?")
    show = counts.rename(index=SHORT_ES)
    st.bar_chart(show.iloc[:, :2], stack=False, horizontal=True, height=280)
    st.dataframe(show, use_container_width=True)
    st.markdown(
        "- **Versión corregida** (por defecto): AICc con los mismos puntos del ajuste, tiempo "
        "medido desde el primer dato y Arrhenius con temperatura de referencia.\n"
        "- **Cálculo de la memoria**: reproduce los scripts originales (N y SST antes del "
        "filtro, C(0) = C₀ aunque la serie parta después, Arrhenius sin T_ref). "
        "Ver `docs/supuestos.md`.")
    st.subheader("Constante de primer orden por sustrato")
    k = res[(res["first_order__k"] > 0) & (res["first_order__k"] < 0.5)]
    fig = go.Figure()
    for g, col in GROUP_COLORS.items():
        d = k[k["Grupo_Residuo"] == g]
        fig.add_trace(go.Box(y=d["first_order__k"], name=GROUP_LABELS[g], marker_color=col,
                             boxpoints="all", jitter=0.4, pointpos=0, text=d.index,
                             hovertemplate="Serie %{text}<br>k = %{y:.3f} d⁻¹<extra></extra>"))
    fig.update_layout(template="plotly_white", yaxis_title="k (d⁻¹)", height=380,
                      showlegend=False, margin=dict(t=10))
    st.plotly_chart(fig, use_container_width=True)
    with st.expander("Tabla completa (159 series)"):
        st.dataframe(res, use_container_width=True)
        st.download_button("Descargar CSV", res.to_csv().encode(), "ajustes_por_serie.csv")

with tab3:
    st.write("Pega o sube una serie con columnas `t` (días), `C` (TOC u OM, %) y, "
             "opcionalmente, `T` (°C) para incluir Arrhenius. Mínimo 5 puntos.")
    ejemplo = ("t,C,T\n0,48.2,25\n4,42.0,58\n9,40.6,62\n15,40.4,55\n21,40.4,45\n"
               "27,40.2,38\n39,39.9,30")
    up = st.file_uploader("CSV", type="csv")
    txt = st.text_area("…o pega aquí", ejemplo, height=180)
    try:
        df_u = pd.read_csv(up) if up else pd.read_csv(pd.io.common.StringIO(txt))
        df_u = df_u.dropna(subset=["t", "C"]).sort_values("t")
        if len(df_u) < 5:
            raise ValueError("se necesitan al menos 5 puntos")
        if (df_u["C"] < 0).any() or (df_u["t"] < 0).any():
            raise ValueError("t y C deben ser no negativos")
        temp = tuple(df_u["T"]) if "T" in df_u and df_u["T"].notna().all() else None
        fits_u = fit_series(tuple(df_u["t"]), tuple(df_u["C"]), temp)
        st.plotly_chart(fit_figure(df_u["t"].to_numpy(), df_u["C"].to_numpy(), fits_u,
                                   "C (%)", "Tus datos"), use_container_width=True)
        st.dataframe(metrics_table(fits_u), hide_index=True, use_container_width=True)
    except Exception as e:  # noqa: BLE001 - se muestra al usuario
        st.error(f"No pude leer los datos: {e}")

with tab_ml:
    st.write("¿Se puede estimar k sin un ensayo cinético? Un bosque aleatorio entrenado con "
             "las 150 series con buen ajuste de primer orden. Ver `docs/modelo_ml.md`.")
    tab_t, model, dom = ml_model()
    c1, c2 = st.columns([1, 2])
    with c1:
        grupo_ml = st.selectbox("Sustrato ", list(GROUP_LABELS), format_func=GROUP_LABELS.get)
        entrada = {
            "Grupo_Residuo": grupo_ml,
            "Init_pH": st.slider("pH inicial", 4.0, 10.0, 7.0, 0.1),
            "Init_CN": st.slider("C/N inicial", 5.0, 80.0, 25.0, 0.5),
            "Init_MC": st.slider("Humedad inicial (%)", 30.0, 85.0, 60.0, 0.5),
            "T_media": st.slider("Temperatura media del ensayo (°C)", 10.0, 75.0, 40.0, 0.5),
        }
    with c2:
        k_pred = ml.predict_k(model, entrada)
        ev = dom.evaluar(entrada)
        st.metric("k predicho (primer orden)", f"{k_pred:.3f} d⁻¹",
                  help="En validación cruzada, 6 de cada 10 predicciones quedan dentro de un "
                       "factor 2 del k medido.")
        if ev["dentro"]:
            st.success("Dentro del dominio de aplicabilidad: hay ensayos parecidos en los datos.")
        else:
            motivo = (f"fuera del rango de: {', '.join(ev['fuera_de_rango'])}"
                      if ev["fuera_de_rango"] else "combinación poco frecuente en los datos")
            st.warning(f"Fuera del dominio de aplicabilidad ({motivo}). La predicción es una "
                       "extrapolación: el error típico se duplica.")
        cv = pd.read_csv(ROOT / "results" / "ml_validacion.csv")
        st.dataframe(cv, hide_index=True, use_container_width=True,
                     column_config={c: st.column_config.NumberColumn(format="%.2f")
                                    for c in ("RMSE_log10", "R2", "dentro_factor_2")})
        st.caption("Conclusión: las condiciones iniciales explican poco de k (R² 0,17 en "
                   "series nuevas y peor que la mediana con un sustrato nuevo). Para diseñar, "
                   "medir la cinética del sustrato real.")

with tab_r:
    st.write("Reactor de compostaje de 32 L con aireación forzada: balances de materia "
             "orgánica, agua, O₂, CO₂ y energía del sólido y del gas (Petric y Selimbašić, "
             "2008). Valores por defecto = experimento 2 del artículo. Ver "
             "`docs/modelo_reactor.md`.")
    c1, c2 = st.columns([1, 3])
    with c1:
        caudal = st.slider("Caudal de aire (m³/h)", 0.02, 0.40, 0.18, 0.01)
        hum = st.slider("Humedad inicial (%)", 40.0, 80.0, 69.1, 0.5)
        masa = st.slider("Masa inicial (kg)", 5.0, 20.0, 12.5, 0.5)
        t0 = st.slider("Temperatura inicial (°C)", 10.0, 40.0, 28.0, 0.5)
        ta = st.slider("Temperatura ambiente (°C)", 5.0, 40.0, 28.0, 0.5)
        corr = st.radio("Factores de corrección de la cinética",
                        ["articulo", "todas"],
                        format_func={"articulo": "Solo temperatura (reproduce el artículo)",
                                     "todas": "T, O₂, humedad y espacio libre"}.get,
                        help="Hallazgo H-1: con los parámetros publicados, solo el factor de "
                             "temperatura reproduce el experimento 2.")
    with c2:
        try:
            sim = simular_reactor(caudal, hum / 100, masa, t0, ta, corr)
        except (ValueError, RuntimeError) as e:
            st.error(str(e))
            st.stop()
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("T máxima", f"{sim['T_C'].max():.1f} °C")
        m2.metric("Horas ≥ 55 °C", f"{int((sim['T_C'] >= 55).sum())} h",
                  help="Referencia de higienización; el criterio normativo no se evalúa aquí.")
        m3.metric("Conversión a 14 d", f"{sim['conversion_pct'].iloc[-1]:.1f} %")
        m4.metric("O₂ mínimo en el gas", f"{sim['O2_pct'].min():.1f} %")
        exp = None
        defecto = (caudal, hum, masa, t0, ta) == (0.18, 69.1, 12.5, 28.0, 28.0)
        if defecto:
            exp = reactor.cargar_experimento(
                ROOT / "data" / "reference" / "petric2008_exp2_digitalizado.csv")
        g1, g2 = st.columns(2)
        paneles = [("T_C", "Temperatura del sustrato (°C)"),
                   ("conversion_pct", "Conversión de materia orgánica (%)"),
                   ("O2_pct", "O₂ en el gas de salida (%, base seca)"),
                   ("humedad_pct", "Humedad del sustrato (%)")]
        for i, (var, titulo) in enumerate(paneles):
            fig = go.Figure()
            fig.add_scatter(x=sim["t_h"], y=sim[var], mode="lines", name="Modelo",
                            line=dict(color="#2a6f97", width=2))
            if exp is not None and var in set(exp["variable"]):
                e = exp[exp["variable"] == var]
                fig.add_scatter(x=e["t_h"], y=e["valor"], mode="markers", name="Experimento 2",
                                marker=dict(color="#c0504d", size=8))
            fig.update_layout(title=titulo, xaxis_title="Tiempo (h)", height=300,
                              margin=dict(l=10, r=10, t=40, b=10), showlegend=exp is not None,
                              legend=dict(orientation="h", y=-0.3))
            (g1 if i % 2 == 0 else g2).plotly_chart(fig, use_container_width=True)
        if defecto:
            st.caption("Puntos: datos digitalizados a mano de las Figs. 3, 4 y 6 del artículo.")
    st.subheader("Efecto del caudal de aire (resto en el experimento 2)")
    st.dataframe(pd.read_csv(ROOT / "results" / "reactor_caudal.csv"), hide_index=True,
                 use_container_width=True,
                 column_config={"caudal_m3_h": st.column_config.NumberColumn("Caudal (m³/h)",
                                                                             format="%.2f"),
                                "T_max_C": st.column_config.NumberColumn("T máx (°C)",
                                                                         format="%.1f"),
                                "t_T_max_h": st.column_config.NumberColumn("Hora de T máx",
                                                                           format="%.0f"),
                                "conversion_final_pct": st.column_config.NumberColumn(
                                    "Conversión (%)", format="%.1f"),
                                "O2_min_pct": st.column_config.NumberColumn("O₂ mín (%)",
                                                                            format="%.1f"),
                                "horas_sobre_55C": st.column_config.NumberColumn(
                                    "Horas ≥ 55 °C", format="%.0f")})
    st.caption("Más aire enfría y seca: baja la temperatura máxima y las horas sobre 55 °C, "
               "pero mantiene el O₂ lejos de condiciones anóxicas. Con poco aire el O₂ cae a ~6 %.")

with tab4:
    st.markdown((ROOT / "docs" / "acerca.md").read_text(encoding="utf-8"))
