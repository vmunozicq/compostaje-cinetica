"""Carga, corrección y segmentación del dataset de compostaje.

Fuente: repositorio público "Compost-Dataset" (licencia MIT),
https://github.com/1191600618/Compost-Dataset, archivo ``Compost dateset.xlsx``
en el commit fijado ``RAW_COMMIT``. Ver ``data/README.md``.

La lógica de segmentación reproduce la de los scripts originales de la memoria
(``procesar_datos.py`` y ``modelacion_masiva.py``).
"""

from __future__ import annotations

import hashlib
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

RAW_COMMIT = "468c14d9bb78fb7804f7e04260b74437e8573caf"
RAW_FILENAME = "Compost dateset.xlsx"
RAW_URL = (
    "https://raw.githubusercontent.com/1191600618/Compost-Dataset/"
    f"{RAW_COMMIT}/Compost%20dateset.xlsx"
)
RAW_SHA256 = "8f32b062344a3de4b85cace67e091a0d3a503547535bdc2f4a3b801b65ecd0f5"

DEFAULT_RAW_PATH = Path(__file__).resolve().parents[2] / "data" / "raw" / RAW_FILENAME

RENAME_MAP = {
    "Waste type": "Waste_Type",
    "Initial pH": "Init_pH",
    "Initial carbon to nitrogen": "Init_CN",
    "Initial moisture content (%)": "Init_MC",
    "Compost time (day)": "Time",
    "总有机碳变化(%)": "TOC_percent",  # carbono orgánico total (%)
    "有机质含量变化(%)": "OM_percent",  # materia orgánica (%)
    "Temperature (℃)": "Temp",
    "Moisture content (%)": "MC",
    "pH": "pH_meas",
}

# Umbral del script original Corregir_Desplazamiento_Columnas.py: una "humedad"
# menor a 20 % es en realidad un pH, lo que delata la rotación de columnas.
ROTATION_MC_THRESHOLD = 20.0

# Límites físico-biológicos por punto (procesar_datos.py). Un NaN no se descarta.
BIO_LIMITS = {
    "Temp": (10, 75),
    "MC": (30, 75),
    "pH_meas": (5.5, 9),
    "C/N": (8, 40),
}

MIN_POINTS = 5
MIN_RANGE = 2.0  # variación mínima de la serie, en puntos porcentuales


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_raw(path: Path = DEFAULT_RAW_PATH) -> Path:
    """Descarga el Excel original desde el commit fijado y verifica su hash."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(RAW_URL, path)
    digest = sha256(path)
    if digest != RAW_SHA256:
        raise ValueError(
            f"El archivo {path} no coincide con la versión fijada "
            f"(sha256 {digest[:12]}… ≠ {RAW_SHA256[:12]}…)."
        )
    return path


def correct_rotated_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Corrige el desplazamiento cíclico pH → Temp → Humedad del archivo original.

    Regla de Vicente (``Corregir_Desplazamiento_Columnas.py``): si la humedad es
    < 20 %, la fila está rotada. En esas filas la columna "pH" trae la
    temperatura, "Temperature" trae la humedad y "Moisture content" trae el pH.
    En la versión fijada del dataset afecta a 1.662 filas (residuos de cocina,
    estiércol y lodos).
    """
    df = df.copy()
    mask = df["MC"] < ROTATION_MC_THRESHOLD
    ph, temp, mc = (df.loc[mask, c].copy() for c in ("pH_meas", "Temp", "MC"))
    df.loc[mask, "Temp"] = ph
    df.loc[mask, "MC"] = temp
    df.loc[mask, "pH_meas"] = mc
    return df


def load_raw(path: Path | str | None = None, correct: bool = True,
             drop_nan_time: bool = True) -> pd.DataFrame:
    """Lee el Excel (fila 0 = grupos, fila 1 = encabezados) y normaliza columnas.

    ``drop_nan_time=False`` reproduce el script de bootstrap de la memoria, que no
    eliminaba la única fila sin tiempo; eso desplaza la numeración de las series
    (ver docs/supuestos.md, D2).
    """
    path = Path(path) if path else download_raw()
    df = pd.read_excel(path, header=1)
    df.columns = [str(c).strip().replace("\n", " ") for c in df.columns]
    df = df.rename(columns=RENAME_MAP)
    if drop_nan_time:
        df = df.dropna(subset=["Time"])
    for c in ("Init_pH", "Init_CN", "Init_MC", "Time", "TOC_percent",
              "OM_percent", "Temp", "MC", "pH_meas", "C/N"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if correct:
        df = correct_rotated_columns(df)
    return df


def assign_series_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Identifica series experimentales.

    1. "Huella" = (tipo de residuo, pH, C/N y humedad iniciales).
    2. Dentro de cada huella, cada vez que el tiempo retrocede empieza una
       réplica nueva ("diente de sierra").
    ID = ``{huella}_{réplica}``, igual que en la memoria (p. ej. ``25_1``).
    """
    df = df.copy()
    group_cols = ["Waste_Type", "Init_pH", "Init_CN", "Init_MC"]
    df[group_cols] = df[group_cols].fillna(-999)
    df["Series_ID_Raw"] = df.groupby(group_cols, sort=True).ngroup() + 1
    parts = []
    for sid, group in df.groupby("Series_ID_Raw", sort=True):
        group = group.sort_index().copy()
        split = group["Time"] < group["Time"].shift(1)
        group["Sub_ID"] = split.cumsum().astype(int) + 1
        group["Unique_ID"] = [f"{sid}_{s}" for s in group["Sub_ID"]]
        parts.append(group)
    return pd.concat(parts)


def apply_bio_limits(df: pd.DataFrame, limits: dict = BIO_LIMITS) -> pd.DataFrame:
    mask = pd.Series(True, index=df.index)
    for col, (lo, hi) in limits.items():
        if col in df.columns:
            mask &= df[col].isna() | df[col].between(lo, hi)
    return df[mask]


def build_validated_dataset(df_raw: pd.DataFrame, min_points: int = MIN_POINTS) -> pd.DataFrame:
    """Equivalente a ``Dataset_Cineticas_Validado_Fisicamente.csv``."""
    df = apply_bio_limits(assign_series_ids(df_raw))
    counts = df["Unique_ID"].value_counts()
    return df[df["Unique_ID"].isin(counts[counts >= min_points].index)].copy()


def waste_group(text: str) -> str:
    t = str(text).lower()
    if "kitchen" in t or "food" in t:
        return "Kitchen Waste"
    if "sludge" in t or "sewage" in t:
        return "Sludge"
    if "manure" in t:
        return "Manure"
    if "green" in t or "yard" in t:
        return "Green Waste"
    return "Other"


@dataclass
class Series:
    uid: str
    waste_type: str
    group: str
    variable: str  # "TOC_percent" u "OM_percent"
    t: np.ndarray
    c: np.ndarray
    temp: np.ndarray
    init_mc: float
    init_cn: float

    @property
    def n(self) -> int:
        return len(self.t)


def extract_series(s: pd.DataFrame, require_temp: bool = True) -> Series | None:
    """Prepara una serie para el ajuste o devuelve None si no pasa los filtros."""
    s = s.sort_values("Time", kind="stable")
    var = "TOC_percent" if s["TOC_percent"].notna().sum() >= MIN_POINTS else "OM_percent"
    cols = ["Time", var, "Temp"] if require_temp else ["Time", var]
    d = s[cols].dropna()
    if len(d) < MIN_POINTS:
        return None
    c = d[var].to_numpy(float)
    if c.max() - c.min() < MIN_RANGE:
        return None
    first = s.iloc[0]
    return Series(
        uid=str(first["Unique_ID"]),
        waste_type=str(first["Waste_Type"]),
        group=waste_group(first["Waste_Type"]),
        variable=var,
        t=d["Time"].to_numpy(float),
        c=c,
        temp=d["Temp"].to_numpy(float) if "Temp" in d else np.full(len(d), np.nan),
        init_mc=float(first["Init_MC"]),
        init_cn=float(first["Init_CN"]),
    )


def iter_series(df_valid: pd.DataFrame, require_temp: bool = True):
    for _, s in df_valid.groupby("Unique_ID", sort=False):
        ser = extract_series(s, require_temp=require_temp)
        if ser is not None:
            yield ser
