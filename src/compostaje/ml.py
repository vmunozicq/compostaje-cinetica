"""Predicción de la constante de primer orden k a partir de las condiciones del ensayo.

Especificación y supuestos en docs/modelo_ml.md (borrador para revisión de Vicente).

Idea: si k se pudiera predecir desde el sustrato y las condiciones iniciales, un
diseñador no necesitaría un ensayo cinético. Se prueba con validación cruzada que
nunca pone la misma serie (ni sus duplicados) en entrenamiento y prueba, y con una
validación "fuera de rango" que deja un sustrato completo fuera del entrenamiento.
Un dominio de aplicabilidad avisa cuando una entrada se aleja de los datos de
entrenamiento, en vez de extrapolar en silencio.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from . import data

NUM_FEATURES = ["Init_pH", "Init_CN", "Init_MC", "T_media"]
CAT_FEATURES = ["Grupo_Residuo"]
TARGET = "log10_k"
# Calidad mínima del ajuste de primer orden para usar su k como dato (docs/modelo_ml.md).
R2_MIN = 0.5
SEED = 0


def build_ml_table(df_raw: pd.DataFrame, res: pd.DataFrame) -> pd.DataFrame:
    """Una fila por serie: condiciones del ensayo, k de primer orden y grupo de CV.

    ``res`` son los ajustes en modo consistente (results/ajustes_por_serie.csv).
    El grupo de validación cruzada es una huella de los datos (t, C): series con
    datos idénticos (duplicados del dataset, supuestos.md L5) quedan en el mismo grupo.
    """
    valid = data.build_validated_dataset(df_raw)
    rows = []
    for s in data.iter_series(valid):
        if s.uid not in res.index:
            continue
        g = valid[valid["Unique_ID"] == s.uid]
        r = res.loc[s.uid]
        rows.append({
            "ID": s.uid, "Grupo_Residuo": s.group,
            "Init_pH": float(g["Init_pH"].iloc[0]), "Init_CN": float(g["Init_CN"].iloc[0]),
            "Init_MC": float(g["Init_MC"].iloc[0]), "T_media": float(np.mean(s.temp)),
            "k": r["first_order__k"], "r2": r["first_order__r2"],
            "huella": hash((tuple(np.round(s.t, 6)), tuple(np.round(s.c, 6)))),
        })
    tab = pd.DataFrame(rows).set_index("ID")
    tab = tab[(tab["r2"] >= R2_MIN) & (tab["k"] > 1e-4)].dropna(subset=NUM_FEATURES)
    tab[TARGET] = np.log10(tab["k"])
    tab["grupo_cv"] = pd.factorize(tab["huella"])[0]
    return tab.drop(columns="huella")


def make_models() -> dict[str, Pipeline]:
    pre = ColumnTransformer([("num", StandardScaler(), NUM_FEATURES),
                             ("cat", OneHotEncoder(handle_unknown="ignore"), CAT_FEATURES)])
    return {
        "Mediana (referencia)": Pipeline([("pre", pre), ("m", DummyRegressor(strategy="median"))]),
        "Ridge": Pipeline([("pre", pre), ("m", RidgeCV(alphas=np.logspace(-2, 3, 20)))]),
        "Bosque aleatorio": Pipeline([("pre", pre), ("m", RandomForestRegressor(
            n_estimators=300, min_samples_leaf=5, random_state=SEED))]),
    }


def _scores(y, yhat) -> dict:
    err = yhat - y
    sst = np.sum((y - y.mean()) ** 2)
    return {"RMSE_log10": float(np.sqrt(np.mean(err ** 2))),
            "R2": float(1 - np.sum(err ** 2) / sst) if sst > 0 else np.nan,
            # Fracción de predicciones dentro de un factor 2 del k ajustado.
            "dentro_factor_2": float(np.mean(np.abs(err) <= np.log10(2)))}


def cross_validate(tab: pd.DataFrame, n_splits: int = 5) -> pd.DataFrame:
    """Dos validaciones por modelo.

    - "Series nuevas": GroupKFold por huella de datos (sustratos conocidos).
    - "Sustrato nuevo": se deja fuera un sustrato completo (fuera de rango).
    """
    x, y = tab[NUM_FEATURES + CAT_FEATURES], tab[TARGET].to_numpy()
    rows = []
    splits = {"Series nuevas": GroupKFold(n_splits=n_splits).split(x, y, tab["grupo_cv"]),
              "Sustrato nuevo": LeaveOneGroupOut().split(x, y, tab["Grupo_Residuo"])}
    for esquema, gen in splits.items():
        folds = list(gen)
        for nombre, model in make_models().items():
            yhat = np.empty_like(y)
            for tr, te in folds:
                model.fit(x.iloc[tr], y[tr])
                yhat[te] = model.predict(x.iloc[te])
            rows.append({"Validación": esquema, "Modelo": nombre, **_scores(y, yhat)})
    return pd.DataFrame(rows)


@dataclass
class DominioAplicabilidad:
    """Distancia media a los k vecinos más cercanos en el espacio estandarizado.

    Una entrada está dentro del dominio si su distancia no supera el percentil 95 de
    las distancias de los propios datos de entrenamiento (excluyéndose a sí mismos), y
    si cada variable numérica cae dentro del rango observado.
    """

    scaler: StandardScaler
    nn: NearestNeighbors
    umbral: float
    rangos: dict
    k: int = 5

    @classmethod
    def fit(cls, tab: pd.DataFrame, k: int = 5) -> DominioAplicabilidad:
        x = tab[NUM_FEATURES].to_numpy(float)
        scaler = StandardScaler().fit(x)
        z = scaler.transform(x)
        nn = NearestNeighbors(n_neighbors=k + 1).fit(z)
        dist, _ = nn.kneighbors(z)
        umbral = float(np.percentile(dist[:, 1:].mean(axis=1), 95))
        rangos = {c: (float(tab[c].min()), float(tab[c].max())) for c in NUM_FEATURES}
        return cls(scaler, nn, umbral, rangos, k)

    def evaluar(self, entrada: dict) -> dict:
        z = self.scaler.transform(np.array([[entrada[c] for c in NUM_FEATURES]], float))
        dist, _ = self.nn.kneighbors(z, n_neighbors=self.k)
        d = float(dist.mean())
        fuera = [c for c in NUM_FEATURES
                 if not self.rangos[c][0] <= entrada[c] <= self.rangos[c][1]]
        return {"distancia": d, "umbral": self.umbral, "fuera_de_rango": fuera,
                "dentro": d <= self.umbral and not fuera}


def fit_final(tab: pd.DataFrame, modelo: str = "Bosque aleatorio"):
    """Entrena el modelo elegido con todos los datos y su dominio de aplicabilidad."""
    model = make_models()[modelo]
    model.fit(tab[NUM_FEATURES + CAT_FEATURES], tab[TARGET].to_numpy())
    return model, DominioAplicabilidad.fit(tab)


def predict_k(model, entrada: dict) -> float:
    """k (d⁻¹) predicho para una entrada con las columnas de NUM_FEATURES y CAT_FEATURES."""
    if not 0 < entrada["Init_MC"] < 100:
        raise ValueError("la humedad inicial debe estar entre 0 y 100 %")
    x = pd.DataFrame([{c: entrada[c] for c in NUM_FEATURES + CAT_FEATURES}])
    return float(10 ** model.predict(x)[0])


def evaluar_dominio(tab: pd.DataFrame, modelo: str = "Bosque aleatorio",
                    n_splits: int = 5) -> pd.DataFrame:
    """¿Sirve el dominio de aplicabilidad para separar predicciones confiables?

    En cada pliegue de "Series nuevas" se entrena el modelo y el dominio solo con el
    pliegue de entrenamiento, y se compara el error de las series de prueba que caen
    dentro y fuera del dominio.
    """
    x, y = tab[NUM_FEATURES + CAT_FEATURES], tab[TARGET].to_numpy()
    rows = []
    for tr, te in GroupKFold(n_splits=n_splits).split(x, y, tab["grupo_cv"]):
        model, dom = fit_final(tab.iloc[tr], modelo)
        yhat = model.predict(x.iloc[te])
        for i, j in enumerate(te):
            ev = dom.evaluar(tab.iloc[j][NUM_FEATURES].to_dict())
            rows.append({"ID": tab.index[j], "Grupo_Residuo": tab.iloc[j]["Grupo_Residuo"],
                         "k_ajustado": 10 ** y[j], "k_predicho": 10 ** yhat[i],
                         "dentro": ev["dentro"], "error_abs_log10": abs(yhat[i] - y[j])})
    return pd.DataFrame(rows).set_index("ID")
