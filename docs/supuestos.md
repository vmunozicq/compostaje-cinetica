# Decisiones, discrepancias y pendientes

La regla del proyecto es que el modelo es mío: el código implementa lo que hice en la memoria. Para eso se compararon tres fuentes: el texto final de la memoria, los scripts de Python que generaron sus tablas y figuras, y esas tablas (copiadas en `data/reference/`). Cuando difieren, el código sigue a los scripts, porque son los que produjeron los números publicados. Cada diferencia queda aquí.

El repositorio tiene dos modos. `modo="memoria"` reproduce los scripts tal cual, y es el que valida `tests/test_validacion.py`. `modo="consistente"`, que es el modo por defecto, corrige D1 y resuelve L1 a L3 (decisión del 2026-11-25).

## Errores encontrados en el cálculo original

**D1. AICc con N y SST de otra versión de la serie.** El script de criterios calculaba SSE = (1 − R²)·SST y usaba N y SST de la serie *antes* del filtro físico-biológico, mientras que el R² venía del ajuste *después* del filtro. En 54 de 159 series N difiere (por ejemplo, 129_1: 6 puntos ajustados, N = 9 en el AICc). El repositorio calcula el AICc con los puntos efectivamente ajustados y conserva el cálculo original en la columna `best_aicc_legacy`:

| Ganador | AICc consistente | Cálculo de la memoria |
| --- | ---: | ---: |
| Primer orden | 80 | 69 (+2)* |
| Gompertz | 59 | 69 |
| Bifásico | 17 | 18 |
| Arrhenius | 1 | 1 |
| Ninguno (todo R² ≤ 0) | 2 | — |

\* El script original rotulaba "1st" a las 2 series sin ningún ajuste válido, por eso la memoria reporta 71.

La conclusión de la memoria se mantiene (no hay modelo universal), pero primer orden gana con más claridad.

**D2. El bootstrap numeró las series distinto.** `analisis_bootstrap.py` no eliminaba la única fila sin tiempo del Excel, lo que corre la numeración. Por eso la serie "25_1" de la Figura 10 (IC 95 % de k: 0,018–0,087 d⁻¹) **no es** la serie 25_1 de la Figura 6. Con la numeración consistente, la 25_1 de la Figura 6 tiene k = 0,34 d⁻¹ (IC 95 % 0,17–0,39). El repositorio usa una sola numeración; `load_raw(drop_nan_time=False)` reproduce la del bootstrap original.

## Afirmaciones de la memoria que los datos no sostienen

Calculadas con las mismas 159 series (idénticas a la tabla original), primero en modo "memoria" y luego en la versión corregida:

- *"Residuos de cocina: mediana de k ≈ 0,12 d⁻¹"*. La mediana es 0,075 d⁻¹ en modo memoria y 0,076 en la versión corregida, con el mismo filtro 0 < k < 0,5. En la versión corregida, la mediana de los tres sustratos queda entre 0,074 y 0,081 d⁻¹.
- *"Los lodos presentan consistentemente λ > 3 d; los residuos de cocina λ → 0"*. En modo memoria, la mediana de λ es 0,94 d en lodos, 1,86 d en cocina y 3,2 d en estiércol. En la versión corregida es 0,85, 1,8 y 0,45 d. El λ alto del estiércol venía de L2.
- *Serie 25_1 como caso de primer orden*. Por AICc gana el bifásico (n = 8), en ambos modos.

**Decisión:** el texto de la memoria no se modifica, porque ya fue entregada. El README, la app y el material de difusión citan los valores corregidos y explican la diferencia.

## Diferencias entre texto y scripts (el repo sigue a los scripts)

**D3. Gompertz.** El Anexo A escribe S(t) = S₀ exp{−exp[(Rm e/S₀)(λ − t) + 1]} como sustrato *remanente*, pero esa función crece de ≈ 0 a S₀. Los scripts la aplican, correctamente, a la degradación acumulada Y = C₀ − C. Así se implementó. El texto del anexo también dice que S → S₀·e⁻¹ para t → ∞, cuando en realidad Y → P.

**D4. Primer orden y bifásico con asíntota.** La memoria escribe S = S₀ e^(−kt), que tiende a 0. Los scripts usan C∞ + (C₀ − C∞) e^(−kt), con C₀ fijo en el primer dato. Así p = 2 (k, C∞) en primer orden y p = 4 (k₁, k₂, f, C∞) en bifásico, como en el cálculo del AICc.

**D5. Cotas.** El Anexo A propone k ∈ [0,001; 0,5] d⁻¹, S₀ ∈ [80, 120] %, A ∈ [0,001; 10] d⁻¹ y Eₐ/R ∈ [1000, 20000] K. Los scripts usan k ∈ [0, 10], C∞ ∈ [0, C₀], Rm ∈ [0, 100], λ ∈ [0, 50], A ∈ [1, 10⁹] y Eₐ ∈ [1, 150] kJ/mol. En el bootstrap: k ≤ 5, Rm ≤ 50, λ ≤ 30. Están centralizadas en `fitting.BOUNDS` y `bootstrap.BOOT_BOUNDS`.

**D6. Filtros de calidad.** El texto habla de eliminar "aumentos repentinos" de materia orgánica y series con degradación neta < 5 %. El script final aplica límites físico-biológicos por punto (T, humedad, pH y C/N), n ≥ 5 y una variación mínima de 2 puntos porcentuales. No hay regla de "aumentos repentinos". El umbral de 5 % solo aparece en un script exploratorio (`modelacion_final_v2.py`, 172 series).

**D7. Temperatura.** El texto menciona interpolar temperaturas faltantes. El script descarta los puntos sin temperatura (para los cuatro modelos, no solo Arrhenius) y solo interpola T(t) entre mediciones dentro del ODE.

**D8. RMSE.** El Anexo B define un RMSE normalizado (%). La figura de la memoria usa RMSE absoluto. El repo entrega ambos (`rmse` y `nrmse_pct`).

**D9. Corrección de columnas del dataset.** En 1.662 filas del archivo original (cocina, estiércol y lodos) las columnas pH, temperatura y humedad están rotadas (p. ej., pH = 29,8). Se aplica mi regla de `Corregir_Desplazamiento_Columnas.py`: si la humedad es < 20 %, se rota. Con eso se recuperan exactamente las 417 series y los 3.362 puntos de la memoria.

## Limitaciones resueltas en la versión corregida

**L1. Arrhenius no era identificable (resuelto).** A y Eₐ se compensan: un A alto con un Eₐ alto da el mismo k. En modo memoria, L-BFGS-B termina en los valores iniciales (A = 10⁴, Eₐ = 30 kJ/mol) en 156 de las 159 series. Por eso los R² de Arrhenius dependen de la versión de SciPy y el test no los compara punto a punto. La versión corregida usa k(T) = k_ref·exp[−Eₐ/R·(1/T − 1/T_ref)], con T_ref igual a la temperatura media medida en la serie. Se ajusta con `least_squares` (TRF) desde Eₐ₀ = 20, 60 y 100 kJ/mol y se conserva el de menor SSE. Mantiene p = 3 y el rango de Eₐ. También se informa el A equivalente.

Con este cambio, el R² mediano sube de 0,85 a 0,96 y Arrhenius gana 10 series, con Eₐ entre 25 y 100 kJ/mol. En 34 series Eₐ queda en la cota inferior (1 kJ/mol): en ellas los datos no muestran efecto de la temperatura.

**L2. Origen del tiempo (resuelto).** C₀ es el primer dato, pero en modo memoria el modelo impone C(0) = C₀ aunque la serie empiece en t > 0. Esto ocurre en 58 de las 159 series; 12 de ellas parten entre los días 5 y 14. La versión corregida mide el tiempo desde el primer dato (t − t₀). Las funciones `predict` siguen recibiendo el tiempo original. Ejemplo: en 461_1, el R² de primer orden sube de 0,90 a 0,98 y k pasa de 0,20 a 0,34 d⁻¹. Gompertz era el más afectado, porque su latencia λ absorbía el desfase.

**L3. Bootstrap de pares con n pequeño (resuelto).** Con 5 a 9 puntos, remuestrear pares produce réplicas con tiempos repetidos y ajustes degenerados. Por ejemplo, el IC de Rm de 454_1 llegaba hasta la cota de 50. Se agregó el bootstrap de residuos: se ajusta una vez, se remuestrean los residuos, se suman a la curva y se reajusta. Es el método por defecto. En 454_1 el IC de Rm queda en 3,9–13,5. Este método supone que el modelo es correcto y que el error es homocedástico (supuesto 5 de `modelo.md`). El de pares se mantiene con `metodo="pares"`.

| Mejor modelo (AICc) | Memoria | Solo D1 | Versión corregida (D1 + L1 + L2) |
| --- | ---: | ---: | ---: |
| Primer orden | 69 (+2) | 80 | 86 |
| Gompertz | 69 | 59 | 40 |
| Bifásico | 18 | 17 | 21 |
| Arrhenius | 1 | 1 | 10 |
| Ninguno | — | 2 | 2 |

## Otras limitaciones

**L4. Semilla.** El bootstrap original no fijaba semilla. Aquí `seed = 42` hace que los IC sean reproducibles, pero no idénticos a los de la memoria.

**L5. Series duplicadas en el dataset.** Algunas series tienen exactamente los mismos ajustes, por ejemplo 11_1 y 135_1, o 43_1 y 430_1. Probablemente son el mismo experimento con dos huellas distintas. **Decisión (2026-11-26):** se mantienen como en la memoria, que las contó por separado.
