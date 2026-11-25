import numpy as np
import pandas as pd

from compostaje import data


def test_correccion_de_columnas_rotadas():
    df = pd.DataFrame({"pH_meas": [55.0, 7.2], "Temp": [60.0, 50.0], "MC": [7.5, 60.0]})
    out = data.correct_rotated_columns(df)
    assert out.loc[0, ["Temp", "MC", "pH_meas"]].tolist() == [55.0, 60.0, 7.5]
    assert out.loc[1, ["Temp", "MC", "pH_meas"]].tolist() == [50.0, 60.0, 7.2]


def test_ids_separan_replicas_por_diente_de_sierra():
    df = pd.DataFrame({"Waste_Type": ["a"] * 6, "Init_pH": 7.0, "Init_CN": 20.0,
                       "Init_MC": 60.0, "Time": [0, 5, 10, 0, 5, 10]})
    ids = data.assign_series_ids(df)["Unique_ID"].tolist()
    assert ids == ["1_1"] * 3 + ["1_2"] * 3


def test_serie_plana_se_descarta():
    s = pd.DataFrame({"Unique_ID": "1_1", "Waste_Type": "food waste", "Init_MC": 60,
                      "Init_CN": 25, "Time": np.arange(6.0), "TOC_percent": 40.0,
                      "OM_percent": np.nan, "Temp": 50.0})
    assert data.extract_series(s) is None
