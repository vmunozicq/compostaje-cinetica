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

from compostaje import bootstrap, data, fitting, pipeline  # noqa: E402
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
def load_results() -> pd.DataFrame:
    return pd.read_csv(ROOT / "results" / "ajustes_por_serie.csv", index_col="ID")


@st.cache_resource(show_spinner="Ajustando modelos…", max_entries=64)
def fit_series(t: tuple, c: tuple, temp: tuple | None):
    # cache_resource: los resultados incluyen funciones de predicción (no serializables).
    return fitting.fit_all(np.array(t), np.array(c), None if temp is None else np.array(temp))


@st.cache_data(show_spinner="Bootstrap en curso…")
def run_bootstrap(t: tuple, c: tuple, n_boot: int):
    out = bootstrap.bootstrap_series(np.array(t), np.array(c), n_boot=n_boot, seed=42)
    return {k: {q: v[q] for q in ("mean", "lo", "hi", "n_ok")} for k, v in out.items()}


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

tab1, tab2, tab3, tab4 = st.tabs(["Explorar una serie", "Resumen de 159 series",
                                  "Ajustar mis datos", "Acerca de"])

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
        n_boot = st.select_slider("Réplicas bootstrap", [100, 250, 500, 1000], value=250)
        do_boot = st.button("Calcular IC 95 % (bootstrap)")
    s = data.extract_series(valid[valid["Unique_ID"] == uid])
    fits = fit_series(tuple(s.t), tuple(s.c), tuple(s.temp))
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
        b = run_bootstrap(tuple(s.t), tuple(s.c), n_boot)
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
    counts = pipeline.winner_counts(res)
    st.subheader("¿Qué modelo gana según AICc?")
    show = counts.rename(index=SHORT_ES)
    st.bar_chart(show.iloc[:, :2], stack=False, horizontal=True, height=280)
    st.dataframe(show, use_container_width=True)
    st.markdown(
        "- **n y SST consistentes**: el AICc se calcula con los mismos puntos usados en el "
        "ajuste (versión por defecto de este repo).\n"
        "- **Cálculo de la memoria**: reproduce el script original, que usaba N y SST de la "
        "serie antes del filtro físico-biológico. Ver `docs/supuestos.md`.")
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

with tab4:
    st.markdown((ROOT / "docs" / "acerca.md").read_text(encoding="utf-8"))
