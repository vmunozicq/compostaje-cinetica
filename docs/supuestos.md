# Decisiones, discrepancias y pendientes

La regla del proyecto es que el modelo es mío: el código implementa lo que hice en la memoria. Para eso se compararon tres fuentes: el texto final de la memoria, los scripts de Python que generaron sus tablas y figuras, y esas tablas (copiadas en `data/reference/`). Cuando difieren, el código sigue a los scripts, porque son los que produjeron los números publicados. Cada diferencia queda aquí.

Las marcadas **TODO(Vicente)** son decisiones abiertas que me toca tomar.

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

Calculadas con las mismas 159 series (idénticas a la tabla original):

- *"Residuos de cocina: mediana de k ≈ 0,12 d⁻¹"*. La mediana es 0,075 d⁻¹ (media 0,104) con el mismo filtro 0 < k < 0,5.
- *"Los lodos presentan consistentemente λ > 3 d; los residuos de cocina λ → 0"*. Mediana de λ: lodos 0,94 d, cocina 1,86 d, estiércol 3,2 d.
- *Serie 25_1 como caso de primer orden*. Por AICc gana el bifásico (n = 8), también con el cálculo original. El primer orden tiene R² = 0,996, pero el bifásico llega a 0,999.

**TODO(Vicente):** decidir si se corrige el texto de esos párrafos en el README y en la presentación, o si se citan los valores nuevos.

## Diferencias entre texto y scripts (el repo sigue a los scripts)

**D3. Gompertz.** El Anexo A escribe S(t) = S₀ exp{−exp[(Rm e/S₀)(λ − t) + 1]} como sustrato *remanente*, pero esa función crece de ≈ 0 a S₀. Los scripts la aplican, correctamente, a la degradación acumulada Y = C₀ − C. Así se implementó. El texto del anexo también dice que S → S₀·e⁻¹ para t → ∞, cuando en realidad Y → P.

**D4. Primer orden y bifásico con asíntota.** La memoria escribe S = S₀ e^(−kt), que tiende a 0. Los scripts usan C∞ + (C₀ − C∞) e^(−kt), con C₀ fijo en el primer dato. Así p = 2 (k, C∞) en primer orden y p = 4 (k₁, k₂, f, C∞) en bifásico, como en el cálculo del AICc.

**D5. Cotas.** El Anexo A propone k ∈ [0,001; 0,5] d⁻¹, S₀ ∈ [80, 120] %, A ∈ [0,001; 10] d⁻¹ y Eₐ/R ∈ [1000, 20000] K. Los scripts usan k ∈ [0, 10], C∞ ∈ [0, C₀], Rm ∈ [0, 100], λ ∈ [0, 50], A ∈ [1, 10⁹] y Eₐ ∈ [1, 150] kJ/mol. En el bootstrap: k ≤ 5, Rm ≤ 50, λ ≤ 30. Están centralizadas en `fitting.BOUNDS` y `bootstrap.BOOT_BOUNDS`.

**D6. Filtros de calidad.** El texto habla de eliminar "aumentos repentinos" de materia orgánica y series con degradación neta < 5 %. El script final aplica límites físico-biológicos por punto (T, humedad, pH y C/N), n ≥ 5 y una variación mínima de 2 puntos porcentuales. No hay regla de "aumentos repentinos". El umbral de 5 % solo aparece en un script exploratorio (`modelacion_final_v2.py`, 172 series).

**D7. Temperatura.** El texto menciona interpolar temperaturas faltantes. El script descarta los puntos sin temperatura (para los cuatro modelos, no solo Arrhenius) y solo interpola T(t) entre mediciones dentro del ODE.

**D8. RMSE.** El Anexo B define un RMSE normalizado (%). La figura de la memoria usa RMSE absoluto. El repo entrega ambos (`rmse` y `nrmse_pct`).

**D9. Corrección de columnas del dataset.** En 1.662 filas del archivo original (cocina, estiércol y lodos) las columnas pH, temperatura y humedad están rotadas (p. ej., pH = 29,8). Se aplica mi regla de `Corregir_Desplazamiento_Columnas.py`: si la humedad es < 20 %, se rota. Con eso se recuperan exactamente las 417 series y los 3.362 puntos de la memoria.

## Limitaciones conocidas y pendientes

**L1. Arrhenius no es identificable.** A y Eₐ se compensan (un A alto con un Eₐ alto da el mismo k). L-BFGS-B con gradiente numérico sobre un ODE converge a mínimos distintos según la versión de SciPy: en 40 de 159 series el R² difiere del original en más de 10⁻³. Por eso el test no compara Arrhenius punto a punto. Gana 1 serie en ambas versiones. **TODO(Vicente):** ¿reparametrizar con una temperatura de referencia, k(T) = k_ref·exp[−Eₐ/R·(1/T − 1/T_ref)], que es la práctica estándar? Cambiaría el modelo, así que no lo hice.

**L2. Origen del tiempo.** C₀ es el primer dato, pero el modelo impone C(0) = C₀ aunque la serie empiece en t > 0 (por ejemplo, 461_1 parte en el día 1). El sesgo es pequeño, pero existe. **TODO(Vicente):** ¿desplazar t → t − t₀?

**L3. Bootstrap de pares con n pequeño.** Con 5 a 9 puntos, remuestrear pares produce réplicas con tiempos repetidos y ajustes degenerados. Algunos IC tocan la cota (454_1: Rm hasta 50). Un bootstrap de residuos sería más estable. **TODO(Vicente):** decidir si se agrega como alternativa.

**L4. Semilla.** El bootstrap original no fijaba semilla. Aquí `seed = 42` hace los IC reproducibles, pero no idénticos a los de la memoria.
