# Ampliación 2: simulador de reactor (calor, agua, O₂ y CO₂)

> **Para revisión de Vicente.** El código está en la rama `reactor` y no se ha integrado a `main`. Cada supuesto o interpretación tiene código (S-R, I) para que puedas aprobarlo o cambiarlo.

## Pregunta de ingeniería

¿Cuánto aire hay que soplar y con qué humedad partir para que el compost se caliente lo suficiente (higienización sobre 55 °C) y degrade la materia orgánica en el menor tiempo? Es la pregunta de la memoria, el tiempo de residencia, respondida con balances de materia y energía en vez de solo un ajuste.

## Fuente

Petric, I., Selimbašić, V. (2008). *Development and validation of mathematical model for aerobic composting process.* Chem. Eng. J. 139(2), 304–317. doi:10.1016/j.cej.2007.08.017. El artículo es de pago: el repositorio cita ecuaciones y constantes, pero no lo redistribuye. Una versión posterior del mismo grupo, con un reactor de 57 L, es Papraćanin (2022), *Kem. Ind.* 71, 185–199 (CC BY 4.0).

¿Por qué no se usan tus 159 series? El Compost Dataset no informa caudal de aire, O₂, volumen del reactor ni temperatura ambiente. Sin esas variables no se puede validar un balance de energía.

## Modelo (números de ecuación del artículo)

Es un reactor de mezcla completa de 32 L con tres fases.

- **Cinética.** dm_OM/dt = −k·m_OM^n (1). k = k_T·k_O₂·k_H₂O·k_FAS (2), con k_T = a[b^(T−20) − c^(T−60)] (4), T en °C. Las correcciones de O₂, humedad y espacio libre son (5), (6) y (7)–(10). Parámetros de la Tabla 5: a = 8,83·10⁻⁵ kg^(1−n) h⁻¹, b = 1,0533, c = 1,2247, n = 2,894.
- **Estequiometría.** Fórmula C₈H₁₉O₅N; Y_O₂ = 0,705, Y_CO₂ = 1,276 e Y_H₂O = 0,360 kg/kg_OM (Tabla 2).
- **O₂ y CO₂ disueltos.** dm_i/dt = R_i − k_La(He_i X_i − p_i) (12)–(14), con Henry según (15)–(16).
- **Agua del material.** dm_w/dt = Y_w·r − k_La(P_s − P_v) (23)–(26).
- **Gases (O₂, CO₂, H₂O, N₂).** Entrada de aire húmedo, transferencia y salida (27)–(31).
- **Energía del gas** (32)–(33) **y del sólido-líquido** (34). Incluye calor de reacción Δh = 1,54·10⁷ J/kg_OM (38), pérdidas por paredes UA = 1900 J h⁻¹ K⁻¹ (36), convección h_c = 10⁵ J h⁻¹ K⁻¹ y entalpía de los gases que salen.

Son 10 variables de estado: masa de materia orgánica, O₂ y CO₂ disueltos, agua, O₂, CO₂, H₂O y N₂ en el gas, temperatura del gas y temperatura del sólido. Se resuelven con `solve_ivp` (LSODA), en `src/compostaje/reactor.py`.

## Supuestos e interpretaciones (para tu aprobación)

- **S-R1.** No se modela el NH₃. Requiere el pH como variable de estado, y el propio artículo lo reproduce mal (su Fig. 7). Su efecto en los balances de O₂ y energía es menor (Y_NH₃ = 0,069).
- **S-R2.** La entalpía del vapor de agua que sale incluye el calor latente (2,44 MJ/kg). El artículo no define h̄_i; sin el latente, el modelo pierde 6,6 kg de agua sin enfriarse, lo que es físicamente imposible.
- **S-R3.** El coeficiente de densidad aparente C de la Ec. (10) no tiene valor en el artículo, solo el rango 0,15–0,4. Se usa 0,25, que deja 16,6 L de gas en el reactor.
- **S-R4.** La temperatura inicial y la ambiente son 28 °C, leídas de la Fig. 3. El artículo no las informa.
- **I-1.** Los c_p de los gases de la Tabla 2 se informan como J kmol⁻¹ K⁻¹, pero sus valores son de J mol⁻¹ K⁻¹. Los del agua, materia orgánica e inorgánica son J kg⁻¹ K⁻¹.
- **I-2.** La Ec. (16), Henry del CO₂, da atm, no Pa: se multiplica por 101 325.
- **I-3.** En (13) y (24) se aplica el signo físico: se consume O₂ y se producen CO₂ y agua.
- **I-4.** El O₂ y el CO₂ se informan en base seca, porque el gas pasa por un condensador antes del analizador.

## Hallazgo H-1: las correcciones de la Ec. 2

Si se aplican las cuatro correcciones de la Ec. (2) con los parámetros publicados, el simulador se queda corto. Llega a 52 °C en vez de 64 °C y convierte el 28 % de la materia orgánica en vez del 39 %. Si se aplica solo k_T, reproduce el experimento y la curva del modelo del propio artículo (máximo cercano a 68 °C en torno a las 28 h).

La lectura más probable es que a, b, c y n se ajustaron con las otras correcciones prácticamente iguales a 1. Esto no se puede confirmar sin el código original. Por eso el simulador usa por defecto `correcciones="articulo"` (solo k_T) y deja `"todas"` como opción. Consecuencia: en la versión por defecto, la humedad y el espacio libre de aire no afectan la velocidad. El caudal de aire sí afecta, a través de la temperatura y del O₂ disponible.

**TODO(Vicente):** confirmar esta elección o, alternativamente, reajustar a, b, c y n con las cuatro correcciones activas. Eso requiere los datos del experimento 1, que el artículo no publica.

## Validación (experimento 2: 12,5 kg, humedad 69,1 %, MO 80,2 % b.s., aire 0,18 m³/h, 14 días)

Los datos medidos se digitalizaron a mano de las Figs. 3, 4 y 6 (`data/reference/petric2008_exp2_digitalizado.csv`). La incertidumbre de lectura es de ±1,5 °C, ±1 punto de conversión y ±0,5 % de O₂.

| Variable | RMSE | Sesgo medio | Error máximo | Modelo del artículo (Tabla 6) |
| --- | ---: | ---: | ---: | --- |
| Temperatura | 5,8 °C | +3,2 °C | 14,0 °C | media 2,3 °C, máx 14,3 °C |
| Conversión de MO | 1,9 puntos | +1,4 | 3,0 | media 1,0, máx 3,8 |
| O₂ de salida | 2,7 % vol | +2,1 | 5,4 | media 0,6, máx 3,3 |

La conversión queda bien reproducida y la temperatura máxima también. Sobre las limitaciones:

- El sistema se enfría más lento que el experimento al final, con unos 4 °C de diferencia. Con S-R4 puede que la temperatura ambiente real haya sido menor.
- La caída de O₂ del primer día y medio es más corta y menos profunda que la medida (mínimo de 16 % contra 12 %).

Criterio de aceptación (`tests/test_reactor.py`): RMSE de conversión < 3 puntos, de temperatura < 7 °C y de O₂ < 3,5 %.

## Resultado de diseño: caudal de aire (`results/reactor_caudal.csv`)

| Caudal (m³/h) | T máx (°C) | Horas sobre 55 °C | Conversión a 14 d | O₂ mínimo |
| ---: | ---: | ---: | ---: | ---: |
| 0,06 | 71 | 144 | 45 % | 6 % |
| 0,10 | 71 | 110 | 43 % | 12 % |
| 0,18 | 69 | 70 | 41 % | 16 % |
| 0,26 | 66 | 47 | 39 % | 18 % |

Más aire enfría y reduce el tiempo sobre 55 °C. Bajo 0,06 m³/h el O₂ se acerca al límite de 5 %, bajo el cual aparecen zonas anaerobias. El artículo encontró su óptimo en 0,10 m³/h, igual a 0,54 L min⁻¹ kg_OM⁻¹. El simulador da la misma dirección, pero con un efecto más suave sobre la conversión.

## Fuera de alcance

Gradientes espaciales (pilas o reactores grandes), NH₃ y pH, y cinética microbiana explícita con biomasa (Monod). El artículo los descarta por ser poco identificables con datos de reactor.
