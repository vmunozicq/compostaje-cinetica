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
