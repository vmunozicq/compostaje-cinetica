"""Figuras estáticas (matplotlib) para el README y el reporte."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .models import MODEL_LABELS  # noqa: E402

INK = "#191919"
INK_2 = "#666666"
GRID = "#e6e6e6"
MODEL_COLORS = {"first_order": "#257de0", "gompertz": "#f47230",
                "biphasic": "#00b894", "arrhenius": "#f2a600"}
MODEL_STYLES = {"first_order": "--", "gompertz": "-", "biphasic": "-.", "arrhenius": ":"}
GROUP_COLORS = {"Kitchen Waste": "#257de0", "Manure": "#f47230", "Sludge": "#00b894"}
GROUP_LABELS = {"Kitchen Waste": "Residuos de cocina", "Manure": "Estiércol", "Sludge": "Lodos"}
SHORT_LABELS = {"1st": "Primer orden", "Gomp": "Gompertz", "Bi": "Bifásico", "Arr": "Arrhenius"}


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)


def plot_series_fits(series, fits, ax=None, title=None):
    ax = ax or plt.subplots(figsize=(7, 4.2))[1]
    t_s = np.linspace(series.t.min(), series.t.max(), 200)
    for m in ("first_order", "biphasic", "gompertz", "arrhenius"):
        r = fits[m]
        if not r.success or not np.isfinite(r.r2):
            continue
        ax.plot(t_s, r.predict(t_s), MODEL_STYLES[m], color=MODEL_COLORS[m], lw=2,
                label=f"{MODEL_LABELS[m]} (R² = {r.r2:.2f})")
    ax.scatter(series.t, series.c, s=40, color=INK, zorder=5, label="Datos experimentales",
               edgecolor="white", linewidth=1.5)
    ax.set_xlabel("Tiempo (d)")
    ax.set_ylabel(f"{'TOC' if series.variable == 'TOC_percent' else 'Materia orgánica'} (%)")
    ax.set_title(title or f"Serie {series.uid}", fontsize=11, loc="left")
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    return ax


def plot_winners(counts, path):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    cols = list(counts.columns[:2])
    idx = [i for i in counts.index if i != "Ninguno"]
    x = np.arange(len(idx))
    w = 0.38
    colors = ["#257de0", "#f47230"]
    for j, c in enumerate(cols):
        vals = counts.loc[idx, c].to_numpy()
        bars = ax.bar(x + (j - 0.5) * w, vals, w - 0.02, color=colors[j], label=c)
        ax.bar_label(bars, fontsize=8, color=INK, padding=2)
    ax.set_xticks(x, [SHORT_LABELS.get(i, i) for i in idx])
    ax.set_ylabel("N.º de series ganadas")
    ax.set_title("Mejor modelo por AICc (159 series)", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_r2_box(res, path):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    models = ["first_order", "biphasic", "gompertz", "arrhenius"]
    data = [res[f"{m}__r2"].dropna().clip(lower=0).to_numpy() for m in models]
    bp = ax.boxplot(data, patch_artist=True, widths=0.5, showfliers=False,
                    medianprops=dict(color=INK, lw=1.5))
    for patch, m in zip(bp["boxes"], models, strict=True):
        patch.set_facecolor(MODEL_COLORS[m])
        patch.set_alpha(0.85)
        patch.set_edgecolor("white")
    ax.set_xticks(range(1, 5), [MODEL_LABELS[m] for m in models])
    ax.set_ylabel("R² (valores negativos = 0)")
    ax.set_title("Bondad de ajuste por modelo", loc="left", fontsize=11)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_k_by_group(res, path, k_max=0.5):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    groups = [g for g in GROUP_COLORS if g in set(res["Grupo_Residuo"])]
    rng = np.random.default_rng(0)
    labels = []
    for i, g in enumerate(groups, start=1):
        k = res.loc[res["Grupo_Residuo"] == g, "first_order__k"].dropna()
        k = k[(k > 0) & (k < k_max)]
        ax.boxplot([k], positions=[i], widths=0.5, showfliers=False,
                   medianprops=dict(color=INK, lw=1.5), boxprops=dict(color=INK_2),
                   whiskerprops=dict(color=INK_2), capprops=dict(color=INK_2))
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(k)), k, s=14, color=GROUP_COLORS[g],
                   alpha=0.8, zorder=3)
        labels.append(f"{GROUP_LABELS[g]}\nmediana {k.median():.3f} · n = {len(k)}")
    ax.set_xticks(range(1, len(groups) + 1), labels)
    ax.set_ylabel("k de primer orden (d⁻¹)")
    ax.set_title(f"Constante cinética por sustrato (0 < k < {k_max})", loc="left", fontsize=11)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_bootstrap(boot_df, path):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    d = boot_df.sort_values("k_mean").reset_index(drop=True)
    y = np.arange(len(d))
    ax.hlines(y, d["k_lo"], d["k_hi"], color="#257de0", lw=2)
    ax.scatter(d["k_mean"], y, color=INK, s=30, zorder=3, edgecolor="white", linewidth=1.5)
    ax.set_yticks(y, [f"{r.ID} · {GROUP_LABELS.get(r.Grupo_Residuo, r.Grupo_Residuo)}"
                      for r in d.itertuples()], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("k de primer orden (d⁻¹, escala log): media e IC 95 % (bootstrap de residuos)")
    ax.set_title("Incertidumbre paramétrica (B = 1000)", loc="left", fontsize=11)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_k_vs(res, xcol, xlabel, path, k_max=0.5):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for g, col in GROUP_COLORS.items():
        d = res[(res["Grupo_Residuo"] == g) & (res["first_order__k"] > 0)
                & (res["first_order__k"] < k_max) & (res[xcol] > 0)]
        ax.scatter(d[xcol], d["first_order__k"], s=22, color=col, label=GROUP_LABELS[g],
                   edgecolor="white", linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("k de primer orden (d⁻¹)")
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_ml_paridad(pred, path):
    """k predicho (validación cruzada) vs. k ajustado, marcando el dominio de aplicabilidad."""
    fig, ax = plt.subplots(figsize=(5.6, 5))
    lims = [pred[["k_ajustado", "k_predicho"]].min().min() / 1.5,
            pred[["k_ajustado", "k_predicho"]].max().max() * 1.5]
    ax.plot(lims, lims, color=INK_2, lw=1)
    ax.fill_between(lims, [v / 2 for v in lims], [v * 2 for v in lims], color=GRID, alpha=0.6,
                    lw=0, label="Dentro de un factor 2")
    for dentro, mk, lab in ((True, "o", "Dentro del dominio"), (False, "X", "Fuera del dominio")):
        d = pred[pred["dentro"] == dentro]
        ax.scatter(d["k_ajustado"], d["k_predicho"], s=30 if dentro else 46, marker=mk,
                   c=[GROUP_COLORS.get(g, INK_2) for g in d["Grupo_Residuo"]],
                   edgecolor="white" if dentro else INK, linewidth=0.8, label=lab, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("k ajustado a los datos (d⁻¹)")
    ax.set_ylabel("k predicho por el bosque aleatorio (d⁻¹)")
    ax.set_title("Predicción de k en series no vistas", loc="left", fontsize=11)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=GROUP_LABELS[g])
               for g, c in GROUP_COLORS.items()]
    ax.legend(handles=ax.get_legend_handles_labels()[0][:1] + handles
              + [plt.Line2D([], [], marker="X", ls="", color=INK, label="Fuera del dominio")],
              frameon=False, fontsize=8, loc="upper left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_reactor_validacion(sims: dict, exp, path):
    """Simulaciones del reactor vs. experimento 2 de Petric y Selimbašić (2008)."""
    paneles = [("T_C", "Temperatura (°C)"), ("conversion_pct", "Conversión de MO (%)"),
               ("O2_pct", "O₂ en el gas de salida (% vol)")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    estilos = [("-", MODEL_COLORS["first_order"]), ("--", MODEL_COLORS["gompertz"])]
    for ax, (var, lab) in zip(axes, paneles, strict=True):
        for (nombre, s), (ls, col) in zip(sims.items(), estilos, strict=False):
            ax.plot(s["t_h"], s[var], ls, color=col, lw=2, label=nombre)
        e = exp[exp["variable"] == var]
        ax.scatter(e["t_h"], e["valor"], s=22, color=INK, zorder=4, label="Experimento",
                   edgecolor="white", linewidth=1)
        ax.set_xlabel("Tiempo (h)")
        ax.set_title(lab, loc="left", fontsize=10)
        ax.set_xticks(range(0, 337, 48))
        _style(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
