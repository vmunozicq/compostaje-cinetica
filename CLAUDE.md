# CLAUDE.md

Instrucciones para Claude Code en este repositorio. Proyecto personal de portafolio de un ingeniero civil químico.

## Reparto de trabajo

- **Vicente escribe el modelo** en `docs/modelo.md`: supuestos, ecuaciones con fuente y caso de validación. Es la especificación.
- **Claude Code escribe** el código, los tests y la interfaz a partir de ese documento.
- Si una ecuación, constante o límite normativo no está en `docs/modelo.md`, **no lo inventes**: pregunta o deja un `TODO(Vicente)` visible. Nunca completes un valor normativo (DS 90, DS 609, DS 46, NCh 1333) de memoria; debe venir citado desde Ley Chile (BCN).
- No copies texto de normas ISO; cita solo número de cláusula.

## Estructura

- `src/compostaje/`: núcleo de cálculo en funciones puras, con unidades en el nombre o en el docstring (`caudal_m3_d`, `conc_mg_L`). Sin dependencias de Streamlit.
- `app/streamlit_app.py`: solo interfaz; llama al núcleo, no recalcula.
- `tests/`: un archivo por módulo, más `test_validacion.py` con el caso de `docs/modelo.md`. La tolerancia del caso de validación se declara en el test con su justificación.
- `data/raw/` es de solo lectura. Todo dato nuevo se registra en `SOURCES.md` (fuente, DOI o URL con commit, licencia, fecha).

## Convenciones

- Python 3.11+. Formato y lint con `ruff`. Tipos en las firmas públicas.
- Comentarios y docstrings en español; nombres de variables en español o inglés técnico, consistentes dentro del módulo.
- Validar entradas físicas (caudales negativos, fracciones fuera de 0 a 1) con `ValueError` y mensaje claro.
- Antes de dar algo por terminado: `ruff check .` y `pytest` en verde.

## Datos

Solo datos públicos o sintéticos. Nada de CMPC, LERV, Claro Vicuña Valenzuela ni COPEC, ni siquiera anonimizado.

## Commits

Mensajes cortos en español, en imperativo ("Agrega carta de control X-barra").
