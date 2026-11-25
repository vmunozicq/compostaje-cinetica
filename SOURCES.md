# Fuentes de datos

| source_id | Fuente | DOI o URL (con commit si es GitHub) | Licencia | Fecha de consulta | Notas |
| --- | --- | --- | --- | --- | --- |
| S01 | "Compost-Dataset" (repositorio de 1191600618 en GitHub), asociado a Zhang, X., Kong, Y., Sun, Z. y Yuan, J. (2025). *Using tree-based machine learning models to predict diverse compost maturity via one-hot encoding*. Waste Management, 205, 114981. | https://github.com/1191600618/Compost-Dataset/blob/468c14d9bb78fb7804f7e04260b74437e8573caf/Compost%20dateset.xlsx | MIT (Copyright (c) 2025 1191600618) | 2026-11-25 | `data/raw/Compost dateset.xlsx`, sin modificar (sha256 `8f32b062…5ecd0f5`). 4.248 muestras compiladas de más de 260 artículos (1993–2025). La corrección de columnas rotadas se aplica en código (`compostaje.data.correct_rotated_columns`). |
| S02 | Muñoz, V. (2025). *Modelación cinética multimodelo del compostaje aeróbico de sólidos a partir de literatura*. Memoria de título, Universidad de Concepción. | — | Del autor (incluida con su permiso) | 2026-11-25 | `data/reference/`: tablas de resultados originales que usa `tests/test_validacion.py`. |

Notas:

- El repositorio de origen publicó en junio de 2026 una versión nueva y estandarizada (`Data_ready.xlsx`, 3.947 registros). Este proyecto usa la versión original de julio de 2025, fijada por commit, que es la que usó la memoria.
- El README del dataset dice que los datos vienen de un "sistema de sensores dedicado". El mismo README también dice que son una compilación de más de 260 artículos, que es como lo describe la memoria. Las referencias primarias por fila no vienen incluidas en el archivo.
- Texto de la licencia MIT de S01: ver [`data/raw/LICENSE-Compost-Dataset`](data/raw/LICENSE-Compost-Dataset).
